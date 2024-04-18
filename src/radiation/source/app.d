import std.math : cos, sin, PI, pow;
import std.stdio;
import std.random : Random, uniform;

import lmr.config;
import lmr.fileutil : ensure_directory_is_present;
import lmr.fluidblock : FluidBlock;
import lmr.fluidfvcell : FluidFVCell;
import lmr.fvcell : FVCell;
import lmr.globalconfig;
import lmr.globaldata;
import lmr.init;
import lmr.sfluidblock;
import lmr.ufluidblock;

import geom.elements.vector3 : Vector3;
import geom.grid.grid;
import nm.number;

import gas.physical_constants : StefanBoltzmann_constant;

/// March across a block from the given cell in desired direction. 
/// This **full** implementation takes many small steps across the domain
/// FIXME: The step size should probably be dependent on cell size...
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

        foreach (neighbour; currentCell.cell_cloud) {
            if (neighbour.is_ghost) {
                // FIXME: Deal with ghost cells!
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
    uint snapshotStart = 0003;
    uint nWrittenSnapshots = snapshotStart + 1;

    alias cfg = GlobalConfig;

    initConfiguration();
    initLocalBlocks();
    initFluidBlocksBasic();
    initFluidBlocksMemoryAllocation();
    initFluidBlocksGlobalCellIDStarts();
    initFluidBlocksZones();
    initFluidBlocksFlowField(snapshotStart);

    auto dirName = snapshotDirectory(nWrittenSnapshots);
    ensure_directory_is_present(dirName);

    auto rng = Random(4);  // Chosen by fair dice roll guaranteed to be random (xkcd.com/221)

    auto decay = 0.1;
    auto firstBlock = localFluidBlocks[0];
    uint angleSamples = 100;

    foreach (cell_id, ref origin; firstBlock.cells) {

        // auto angle = uniform(0.0f, 2*PI, rng);
        number fullEmission = StefanBoltzmann_constant * origin.fs.gas.T ^^ 4;

        foreach (a; 0 .. angleSamples) {
            // Loop through angles for now
            auto angle = 2*PI * a / angleSamples;
            auto direction = Vector3([cos(angle), sin(angle), 0]);
            direction.normalize(); // Shouldn't be needed with random angle

            auto rayStrength = fullEmission / angleSamples;
            origin.fs.Qrad -= rayStrength;  // Remove the energy of the ray emitted

            size_t[] rayCells;
            number[] rayLengths;
            marching_full(cell_id, firstBlock, direction, rayCells, rayLengths);
            writeln("Stepped over ", rayCells.length, " cells");

            number heating;
            foreach (i; 0 .. rayCells.length) {
                heating = rayStrength * (1 - (1 - decay) ^^ rayLengths[i]);
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

