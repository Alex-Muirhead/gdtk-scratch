import std.stdio;

import lmr.config;
import lmr.fileutil : ensure_directory_is_present;
import lmr.fvcell;
import lmr.globalconfig;
import lmr.globaldata;
import lmr.init;
import lmr.ufluidblock;
import lmr.sfluidblock;

import geom.elements.vector3;
import geom.grid.grid;

import gas.physical_constants : StefanBoltzmann_constant;

void main() {
    writeln("Initial workings on a standalone radiation post-processing code.");

    // Hard-code how many snapshots we're working with
    auto snapshotStart = 0003;

    alias cfg = GlobalConfig;

    initConfiguration();
    initLocalBlocks();
    initFluidBlocksBasic();
    initFluidBlocksMemoryAllocation();
    initFluidBlocksGlobalCellIDStarts();
    initFluidBlocksZones();
    initFluidBlocksFlowField(snapshotStart);

    int nWrittenSnapshots = 4;

    Grid[] blockGrids = [];
    foreach (blk; localFluidBlocks) {
        foreach (cell; blk.cells) {
            cell.fs.Qrad = StefanBoltzmann_constant * cell.fs.gas.T ^^ 4;
        }

        final switch (blk.grid_type) {
            case Grid_t.structured_grid:
                SFluidBlock sblk = cast(SFluidBlock) blk;
                blockGrids ~= sblk.grid;
                break;
            case Grid_t.unstructured_grid:
                UFluidBlock ublk = cast(UFluidBlock) blk;
                blockGrids ~= ublk.grid;
                break;
        }
    }
    writeln("Grid? ", blockGrids);

    auto dirName = snapshotDirectory(nWrittenSnapshots);
    ensure_directory_is_present(dirName);

    size_t cell_id = 0;

    auto firstBlock = localFluidBlocks[0];
    auto cell = firstBlock.cells[cell_id];
    auto decay = 0.01; // NOTE: Should be proportional to the step size

    Vector3 rayPoint = cell.pos[0];
    auto rayStrength = cell.fs.Qrad;

    // Structure doesn't move, so we can take first position
    writeln("Cell 0 at pos: ", rayPoint);

    auto direction = Vector3([1, 1, 0]);
    direction.normalize();
    auto stepLength = 0.001;

    for (auto i = 0; i < 1000; i++) {
        rayPoint += direction * stepLength;
        auto heating = decay * rayStrength;
        rayStrength -= heating;

        bool contained = blockGrids[0].point_is_inside_cell(rayPoint, cell_id);
        if (!contained) {
            size_t neighbour_id;
            inner: foreach (neighbour; cell.cell_cloud) {
                neighbour_id = neighbour.id;
                if (neighbour_id > 1_000_000_000) {
                    // HACK: Dealing with weird ghost-cell ids
                    continue inner;
                }
                if (blockGrids[0].point_is_inside_cell(rayPoint, neighbour_id)) {
                    contained = true;
                    break inner;
                }
            }
            if (contained) {
                // Successfully found next point
                cell_id = neighbour_id;
                cell = firstBlock.cells[cell_id];
            } else {
                writeln("Error! Couldn't find the thing");
                break;
            }
        }
        writeln("Moved to point: ", rayPoint, " within cell ", cell_id);

        // This has to be at the end, since we have to find *where* to dump the 
        // energy before ending the routine. 
        if (rayStrength < 1E-3) { // Some minimum energy threshold
            break;
        }
    }

    // The method `geom.grid.usgrid.UnstructuredGrid.get_list_of_boundary_cells`
    // looks at what is inside a cell vs outside
    // There is also the method `geom.grid.grid.Grid.point_is_inside_cell`
    // How do we get access to the grid though?
    
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

