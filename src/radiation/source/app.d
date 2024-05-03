import std.algorithm.iteration;
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
    size_t cellID, FluidBlock block, Vector3 direction, 
    ref size_t[] crossed, ref number[] lengths
) {
    FluidFVCell currentCell = block.cells[cellID];
    FVInterface outgoing;
    Vector3 rayCoord = currentCell.pos[0];  // INFO: Grid is stationary, index is time
    Grid currentGrid = get_grid(block);

    
    foreach (bigIndex; 1..10) {
        // INFO: Grid is stationary, index is time
        // FIXME: We should probably operate on the grid directly when integrating into lmr
        Vector3[] vertices = currentCell.vtx.map!(item => item.pos[0]).array;

        size_t nvtx; size_t[8] vtx_id;
        currentGrid.copy_vtx_id_list_for_cell(vtx_id, nvtx, cellID);

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

        dim: switch (currentGrid.dimensions) {
        case 1: throw new Exception("cell search not implemented for 1D grids");
        case 2:
            shape: switch (nvtx) {
            // case 3:
            //     // 3 sets of `on_left_of_xy_line`
            //     inside_cell = inside_xy_triangle(vertices[vtx_id[0]], vertices[vtx_id[1]],
            //                                      vertices[vtx_id[2]], p);
            //     break;
            case 4:
                // 4 sets of `on_left_of_xy_line`
                Vector3 p0 = currentGrid.vertices[vtx_id[0]];
                Vector3 p1 = currentGrid.vertices[vtx_id[1]];
                Vector3 p2 = currentGrid.vertices[vtx_id[2]];
                Vector3 p3 = currentGrid.vertices[vtx_id[3]];
                Vector3 q = rayCoord;

                Vector3 p0diff = q - p0;
                Vector3 p1diff = q - p1;
                Vector3 p2diff = q - p2;
                Vector3 p3diff = q - p3;

                Vector3 p01 = p1 - p0;
                Vector3 p12 = p2 - p1;
                Vector3 p23 = p3 - p2;
                Vector3 p30 = p0 - p3;

                Vector3 d = direction;

                // Coordinate along edge
                number s01 = wedge2D(p0diff, d) / wedge2D(p01, d);
                number s12 = wedge2D(p1diff, d) / wedge2D(p12, d);
                number s23 = wedge2D(p2diff, d) / wedge2D(p23, d);
                number s30 = wedge2D(p3diff, d) / wedge2D(p30, d);

                // Length that the ray has stepped
                // Ensure is positive, to make us move fowards
                number t01 = wedge2D(p0diff, p01) / wedge2D(p01, d);
                number t12 = wedge2D(p1diff, p12) / wedge2D(p12, d);
                number t23 = wedge2D(p2diff, p23) / wedge2D(p23, d);
                number t30 = wedge2D(p3diff, p30) / wedge2D(p30, d);

                size_t boundary_id;
                number step_length;
                if (s01 >= 0 && s01 < 1 && t01 > 0) { 
                    boundary_id = 0; 
                    step_length = t01;
                }
                if (s12 >= 0 && s12 < 1 && t12 > 0) { 
                    boundary_id = 1; 
                    step_length = t12;
                }
                if (s23 >= 0 && s23 < 1 && t23 > 0) { 
                    boundary_id = 2; 
                    step_length = t23;
                }
                if (s30 >= 0 && s30 < 1 && t30 > 0) { 
                    boundary_id = 3; 
                    step_length = t30;
                }

                rayCoord += direction * step_length;

                break shape;
            default:
                throw new Exception("invalid cell type in 2D");
            } // end switch (vtx_id.length)
            break dim;
        case 3:
            switch (nvtx) {
            // case 4:
            //     // 4 sets of `tetrahedron_volume(...) < 0`
            //     inside_cell = inside_tetrahedron(vertices[vtx_id[0]], vertices[vtx_id[1]],
            //                                      vertices[vtx_id[2]], vertices[vtx_id[3]], p);
            //     break;
            // case 8:
            //     // 6 sets of `tetragonal_dypyramid_volume(...) < 0`
            //     inside_cell = inside_hexahedron(vertices[vtx_id[0]], vertices[vtx_id[1]],
            //                                     vertices[vtx_id[2]], vertices[vtx_id[3]],
            //                                     vertices[vtx_id[4]], vertices[vtx_id[5]],
            //                                     vertices[vtx_id[6]], vertices[vtx_id[7]], p);
            //     break;
            // case 5:
            //     inside_cell = inside_pyramid(vertices[vtx_id[0]], vertices[vtx_id[1]],
            //                                  vertices[vtx_id[2]], vertices[vtx_id[3]],
            //                                  vertices[vtx_id[4]], p);
            //     break;
            // case 6:
            //     inside_cell = inside_wedge(vertices[vtx_id[0]], vertices[vtx_id[1]],
            //                                vertices[vtx_id[2]], vertices[vtx_id[3]],
            //                                vertices[vtx_id[4]], vertices[vtx_id[5]], p);
            //     break;
            default:
                throw new Exception("invalid cell type in 3D");
            } // end switch (vtx_id.length)
            break dim;
        default: throw new Exception("invalid number of dimensions");
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
    ref size_t[] crossed, ref number[] lengths) 
{
    bool isWithinBlock = true;
    bool isContained;
    FluidFVCell currentCell = block.cells[cellID];
    Vector3 rayCoord = currentCell.pos[0];  // INFO: Grid is stationary, index is time
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
                    hit = currentCell.iface[i-1]; // 0th cell is self
                    break marchLoop;
                }
            } else {
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

    version(mpi_parallel) { MPI_Barrier(MPI_COMM_WORLD); }

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

    auto rng = Random(4);  // Chosen by fair dice roll guaranteed to be random (xkcd.com/221)

    auto decay = 0.5;
    auto firstBlock = localFluidBlocks[0];
    uint angleSamples = 100;

    number baseEmission = StefanBoltzmann_constant * 280 ^^ 4;

    foreach (cell_id, ref origin; firstBlock.cells) {

        // auto angle = uniform(0.0f, 2*PI, rng);
        number fullEmission = origin.volume[0] * StefanBoltzmann_constant * origin.fs.gas.T ^^ 4;
        fullEmission -= baseEmission; // Add back baseline background radiation
        origin.fs.Qrad -= fullEmission;  // Remove the energy of the ray emitted

        foreach (a; 0 .. angleSamples) {
            // Loop through angles for now
            auto angle = 2*PI * a / angleSamples;
            auto direction = Vector3([cos(angle), sin(angle), 0]);
            direction.normalize(); // Shouldn't be needed with random angle

            auto rayStrength = fullEmission / angleSamples;

            size_t[] rayCells;
            number[] rayLengths;
            FVInterface hit = marching_full(cell_id, firstBlock, direction, rayCells, rayLengths);
            writeln("Interface: ", hit);

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
                heating = rayStrength * (1 - (1-decay) ^^ rayLengths[i]);
                rayStrength -= heating;
                firstBlock.cells[rayCells[i]].fs.Qrad += heating;
            }
        }
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

