import std.stdio;

import fileutil; // Also from lmr
import globalconfig; // Also from lmr
import globaldata; // Also from lmr
import lmr.fvcell;
import init; // Also from lmr
import lmrconfig;

void main() {
    writeln("Initial workings on a standalone radiation post-processing code.");

    auto snapshotStart = 0003;

    alias cfg = GlobalConfig;

    initConfiguration();
    initLocalBlocks();
    initFluidBlocksBasic();
    initFluidBlocksMemoryAllocation();
    initFluidBlocksGlobalCellIDStarts();
    initFluidBlocksZones();
    initFluidBlocksFlowField(snapshotStart);

    foreach (blk; localFluidBlocks) {
        foreach (cell; blk.cells) {
            writeln("Cell id: ", cell.id);
            writeln("Cell fs: ", cell.fs);
        }
    }

    int nWrittenSnapshots = 4;

    auto dirName = snapshotDirectory(nWrittenSnapshots);
    ensure_directory_is_present(dirName);

    // Will currently break, since we haven't constructed the mCIO!
    // Huh... broke because the directly didn't exist
    // ...And with that directory fix, we're fine!
    // So where is it getting the mCIO data from?
    // Ah! It must be from when we call `initFluidBlocksFlowField`

    foreach (blk; localFluidBlocks) {
        auto fileName = fluidFilename(nWrittenSnapshots, blk.id);
        FVCell[] cells;
        cells.length = blk.cells.length;
        foreach (i, ref c; cells) c = blk.cells[i];
        fluidBlkIO.writeVariablesToFile(fileName, cells);
    }

    // So is it possible to add a new variable to this list?
}

