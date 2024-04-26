import std.math : cos, sin, PI, pow;
import std.stdio;
import std.random : Random, uniform;

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

import geom.elements.vector3 : Vector3;
import geom.grid.grid;
import nm.number;

import gas.physical_constants : StefanBoltzmann_constant;

struct Ray {
    Vector3 position;
    Vector3 direction;
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
 */
void marching_full(
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
                // FIXME: Deal with ghost cells!
                FVInterface inter = currentCell.iface[i-1]; // 0th cell is self
                if (!inter.is_on_boundary) {
                    // throw new Exception("Ray hit a ghost cell not attached to a boundary");
                    continue;
                }

                BoundaryCondition boundary = block.bc[inter.bc_id];
                auto mygce = cast (GhostCellFullFaceCopy) boundary.preReconAction[0];
                if (mygce) {
                    // We have a full face copy, does this only exist for structured?
                    // I think we use `GhostCellMappedCopy` for unstructured...
                    // QUESTION: Are there other effects that mean we walk across the boundary?

                    // NOTE: We should always intersect the 0th (first) layer of the ghost cells
                    //       in the boundary. Since ghost cells are added walking away from the
                    //       boundary, this means we can multiply the `i_bndry` by the number
                    //       of ghost cell layers to get the correct position of the ghost cells
                    //       and consequently the mapped cell. 
                    size_t boundaryCellID = inter.i_bndry * block.n_ghost_cell_layers;
                    FluidFVCell mappedCell = mygce.mapped_cells[boundaryCellID];
                    writeln("Cell ", currentCell.id, " maps to ", mappedCell.id);
                }
                continue;
            }
            isContained = currentGrid.point_is_inside_cell(rayCoord, neighbour.id);
            if (isContained) {
                currentCell = neighbour;
                continue marchLoop;
            }
        }

        if (!isContained) {
            // We haven't found the neighbour
            break marchLoop;
        }
    }
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

void main() {
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
            marching_full(cell_id, firstBlock, direction, rayCells, rayLengths);

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

