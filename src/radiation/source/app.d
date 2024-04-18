import std.math : cos, sin, PI;
import std.stdio;
import std.random : Random, uniform;

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
            // Set everything to zero to be sure
            cell.fs.Qrad = 0.0;
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

    auto rng = Random(4);  // Chosen by fair dice roll guaranteed to be random (xkcd.com/221)

    auto decay = 0.01; // NOTE: Should be proportional to the step size
    auto stepLength = 0.001;

    auto firstBlock = localFluidBlocks[0];

    foreach (cell_id, ref origin; firstBlock.cells) {

        // auto angle = uniform(0.0f, 2*PI, rng);

        for (auto a = 0.0; a < 1; a += 0.01) {
            // Loop through angles for now
            auto angle = 2*PI * a;
            auto direction = Vector3([cos(angle), sin(angle), 0]);
            direction.normalize(); // Shouldn't be needed with random angle

            writeln("Moving at angle: ", angle);

            auto cell = origin;
            Vector3 rayPoint = cell.pos[0]; // Grid isn't moving, take the first position

            auto rayStrength = StefanBoltzmann_constant * cell.fs.gas.T ^^ 4;
            cell.fs.Qrad -= rayStrength;  // Remove the energy of the ray emitted

            steps: foreach (i; 0 .. 1000) {
                rayPoint += direction * stepLength;
                auto heating = decay * rayStrength;
                rayStrength -= heating;

                bool contained = blockGrids[0].point_is_inside_cell(rayPoint, cell_id);
                if (!contained) {
                    size_t neighbour_id;
                    inner: foreach (neighbour; cell.cell_cloud) {
                        neighbour_id = neighbour.id;
                        Vector3[][] vertices = [];
                        foreach (j, ref vertex; neighbour.vtx) {
                            vertices ~= vertex.pos;
                        }
                        if (neighbour.is_ghost) {
                            // FIXME: Deal with ghost cells at the boundaries
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
                        // Could not find next cell
                        // writeln("Error! Couldn't find the thing");
                        break steps;
                    }
                }
                cell.fs.Qrad += heating * 10;

                // This has to be at the end, since we have to find *where* to dump the 
                // energy before ending the routine. 
                if (rayStrength < 1E-3) { // Some minimum energy threshold
                    break steps;
                }
            }
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

