import std.algorithm.mutation;
import std.array;
import std.math : cos, sin, PI, pow;
import std.stdio;
import std.random : Random, uniform;
import std.getopt;

import lmr.config;
import lmr.bc.boundary_condition : BoundaryCondition;
import lmr.bc.ghost_cell_effect.full_face_copy : GhostCellFullFaceCopy;
import lmr.fileutil : ensure_directory_is_present;
import lmr.fluidblock : FluidBlock;
import lmr.fluidfvcell : FluidFVCell;
import lmr.fvcell : FVCell;
import lmr.fvinterface : FVInterface;
import lmr.globalconfig;
import lmr.globaldata;
import lmr.init;
import lmr.sfluidblock;
import lmr.ufluidblock;

import geom.elements.vector3 : Vector3, wedge2D;
import geom.grid.grid;
import nm.number;

import gas.physical_constants : StefanBoltzmann_constant;

struct Ray {
    Vector3 position;
    Vector3 direction;
}

FVInterface marching_efficient(
    size_t cellID, FluidBlock block, Vector3 rayTangent,
    ref size_t[] crossed, ref number[] lengths
) {
    // Placeholder for if we need to do moving grids
    size_t gtl = 0;

    FluidFVCell currentCell = block.cells[cellID];
    Vector3 rayCoord = currentCell.pos[gtl];
    Grid currentGrid = get_grid(block);
    bool isWithinBlock = true;

    // Working variables
    FVInterface outgoing;
    Vector3 vertex_i;
    Vector3 vertex_j;

    // NOTE: Push all this down into the grid level
    //       Potentially with the following signature

    /**
     * Params:
     *   line = A description of the ray position & orientation
     *   indx = The cell index to check for intersections at (note, is this enough info?)
     *   crossings = The interfaces which the line crosses through. The two (2) solutions correspond
     *               to the nonpositive (<=0) and positive (>0) solutions.
     */
    // bool compute_line_intersections(ref Ray line, ref size_t indx, out FVInterface[2] crossings)

    rayTangent.normalize();

    while (isWithinBlock) {
        ndim: switch (currentGrid.dimensions) {
        case 1:
            throw new Exception("cell search not implemented for 1D grids");
        case 2:
            size_t iface_id;
            number step_length;

            faces: foreach (n, iface; currentCell.iface) {
                vertex_i = iface.vtx[0].pos[gtl];
                vertex_j = iface.vtx[1].pos[gtl];
                if (currentCell.outsign[n] == -1) {
                    swap(vertex_i, vertex_j);
                } 

                Vector3 toVertex = rayCoord - vertex_i;
                // NOTE: Could potentially use iface.t1 for this?
                Vector3 faceTangent = vertex_j - vertex_i;

                number alignment = wedge2D(faceTangent, rayTangent);
                number s = wedge2D(toVertex, faceTangent) / alignment;
                number t = wedge2D(toVertex, rayTangent) / alignment;

                // There should only be a single solution if the cell is convex
                if (t >= 0 && t < 1 && s > number.epsilon) {
                    iface_id = n;
                    step_length = s;
                    break faces;
                }
            }

            rayCoord += rayTangent * step_length;
            outgoing = currentCell.iface[iface_id];

            lengths ~= step_length;
            crossed ~= currentCell.id;

            // WARN: The `currentCell` can be null if crossing a boundary
            //       with no ghost cell on the other side
            currentCell = (currentCell.outsign[iface_id] == +1)
                ? outgoing.right_cell
                : outgoing.left_cell;

            if (outgoing.is_on_boundary) {
                isWithinBlock = false;
            }

            break ndim;
        case 3:
            throw new Exception("cell search not implemented for 3D grids");
        default:
            throw new Exception("invalid number of dimensions");
        } // end switch (dimensions)
    }

    return outgoing;
}

// FIXME: The step size should probably be dependent on cell size...
// NOTE: What do we want to return here?
//       I'm thinking that we could return the FVInterface or BoundaryCondition
//       that the ray hits, and then make a decision after that.
//       This allows us to deal with different conditions (i.e. periodic boundaries, etc)
//       and also work out what the next block is to restart the routine.

/** March across a block from the given cell in desired direction.
 * This **full** implementation takes many small steps across the domain
 *
 * Params:
 *   cellID = The ID of the starting cell (must be within the block)
 *   block = The block to march the ray across
 *   direction = The direction of the ray as a Vector3
 *   crossed = Which cells within the block the ray crosses
 *   lengths = What distance across each respective cell the ray marched
 * Returns:
 *   The interface which was intersected
 */
FVInterface marching_full(
    size_t cellID, FluidBlock block, Vector3 direction,
    ref size_t[] crossed, ref number[] lengths) {
    bool isWithinBlock = true;
    bool isContained;
    FluidFVCell currentCell = block.cells[cellID];
    Vector3 rayCoord = currentCell.pos[0]; // INFO: Grid is stationary, index is time
    Grid currentGrid = get_grid(block);

    number stepSize = 1E-03; // FIXME: This should be proportional to something...
    int consecutiveSteps = 1;
    FVInterface hit;

    // NOTE: Do we want to record the first step?

    marchLoop: while (isWithinBlock) {
        rayCoord += stepSize * direction;
        isContained = currentGrid.point_is_inside_cell(rayCoord, currentCell.id);

        if (isContained) {
            consecutiveSteps++;
            continue;
        }

        // Record current steps
        crossed ~= currentCell.id;
        lengths ~= stepSize * consecutiveSteps;
        consecutiveSteps = 1;

        // NOTE: Could shortcut by checking the boundaries directly from the cell
        //       vertices, rather than bringing all neighbours into memory
        foreach (i, neighbour; currentCell.cell_cloud) {
            if (neighbour.is_ghost) {
                // Need a different way to check if we're inside the ghost cell
                size_t nvtx = neighbour.vtx.length;
                writeln("Ghost cell has ", nvtx, " vertices.");
                if (true) {
                    hit = currentCell.iface[i - 1]; // 0th cell is self
                    break marchLoop;
                }
            }
            else {
                isContained = currentGrid.point_is_inside_cell(rayCoord, neighbour.id);
                if (isContained) {
                    currentCell = neighbour;
                    continue marchLoop;
                }
            }
        }

        if (!isContained) {
            // We haven't found the neighbour
            break marchLoop;
        }
    }
    return hit;
}

Grid get_grid(FluidBlock block) {
    final switch (block.grid_type) {
    case Grid_t.structured_grid:
        SFluidBlock sblk = cast(SFluidBlock) block;
        return sblk.grid;
    case Grid_t.unstructured_grid:
        UFluidBlock ublk = cast(UFluidBlock) block;
        return ublk.grid;
    }
}

void main(string[] args) {

    string workingDir = ".";
    getopt(args,
        std.getopt.config.stopOnFirstNonOption,
        "d|dir", &workingDir
    );

    writeln("Initial workings on a standalone radiation post-processing code.");

    // Hard-code how many snapshots we're working with
    uint snapshotStart = 0004;
    uint nWrittenSnapshots = snapshotStart + 1;

    alias cfg = GlobalConfig;

    initConfiguration();
    initLocalBlocks();
    initFluidBlocksBasic();
    initFluidBlocksMemoryAllocation();
    initFluidBlocksGlobalCellIDStarts();
    initFluidBlocksZones();
    initFluidBlocksFlowField(snapshotStart);

    version (mpi_parallel) {
        MPI_Barrier(MPI_COMM_WORLD);
    }

    initFullFaceDataExchange();
    initMappedCellDataExchange();
    initGhostCellGeometry();

    // Set everything to zero initially
    // FIXME: This is because of some weird buffer thing in loading
    //        data during `lmr snapshot2vtk`, it takes the same
    //        values from the previous variable (temperature here)
    foreach (blk; localFluidBlocks) {
        foreach (cell; blk.cells) {
            cell.fs.Qrad = 0.0;
        }
    }

    auto dirName = snapshotDirectory(nWrittenSnapshots);
    ensure_directory_is_present(dirName);

    auto rng = Random(4); // Chosen by fair dice roll guaranteed to be random (xkcd.com/221)

    auto decay = 0.5;
    auto firstBlock = localFluidBlocks[0];
    uint angleSamples = 8;

    number baseEmission = StefanBoltzmann_constant * 280 ^^ 4;

    outer: foreach (cell_id, ref origin; firstBlock.cells) {

        // auto angle = uniform(0.0f, 2*PI, rng);
        number fullEmission = origin.volume[0] * StefanBoltzmann_constant * origin.fs.gas.T ^^ 4;
        fullEmission -= baseEmission; // Add back baseline background radiation
        origin.fs.Qrad -= fullEmission; // Remove the energy of the ray emitted

        foreach (a; 0 .. angleSamples) {
            // Loop through angles for now
            auto angle = 2 * PI * a / angleSamples;
            auto direction = Vector3([cos(angle), sin(angle), 0]);
            direction.normalize(); // Shouldn't be needed with random angle

            auto rayStrength = fullEmission / angleSamples;

            size_t[] rayCells;
            number[] rayLengths;
            // FVInterface hit = marching_full(cell_id, firstBlock, direction, rayCells, rayLengths);
            FVInterface hit = marching_efficient(cell_id, firstBlock, direction, rayCells, rayLengths);
            writeln("Interface: ", hit);
            continue;

            //     if (!inter.is_on_boundary) {
            //         writeln("Neighbouring cell: ", neighbour);
            //         throw new Exception("Ray hit a ghost cell not attached to a boundary");
            //         continue;
            //     }

            //     BoundaryCondition boundary = block.bc[inter.bc_id];
            //     auto mygce = cast (GhostCellFullFaceCopy) boundary.preReconAction[0];
            //     if (mygce) {
            //         // We have a full face copy, does this only exist for structured?
            //         // I think we use `GhostCellMappedCopy` for unstructured...
            //         // QUESTION: Are there other effects that mean we walk across the boundary?

            //         // NOTE: We should always intersect the 0th (first) layer of the ghost cells
            //         //       in the boundary. Since ghost cells are added walking away from the
            //         //       boundary, this means we can multiply the `i_bndry` by the number
            //         //       of ghost cell layers to get the correct position of the ghost cells
            //         //       and consequently the mapped cell.
            //         size_t boundaryCellID = inter.i_bndry * block.n_ghost_cell_layers;
            //         FluidFVCell mappedCell = mygce.mapped_cells[boundaryCellID];
            //         writeln("Cell ", currentCell.id, " maps to ", mappedCell.id);
            //     }
            //     continue;
            // }
            number heating;
            foreach (i; 0 .. rayCells.length) {
                heating = rayStrength * (1 - (1 - decay) ^^ rayLengths[i]);
                rayStrength -= heating;
                firstBlock.cells[rayCells[i]].fs.Qrad += heating;
            }
        }
        break outer;
    }

    foreach (blk; localFluidBlocks) {
        auto fileName = fluidFilename(nWrittenSnapshots, blk.id);
        FVCell[] cells;
        cells.length = blk.cells.length;
        //       +-- Enumeration count
        //       v  vvv-- Makes c a pointer, so we can assign back to cells
        foreach (i, ref c; cells) {
            c = blk.cells[i];
        }
        fluidBlkIO.writeVariablesToFile(fileName, cells);
    }
}
