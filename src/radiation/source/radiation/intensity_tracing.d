import std.algorithm.mutation : swap;
import std.algorithm : min, max;
import std.array;
import std.format;
import std.math : cos, sin, PI, pow, fabs;
import std.stdio;
import std.getopt;

import geom.elements.nomenclature : Face, opposite_face;
import geom.elements.vector3 : Vector3, distance_between;
import gas.physical_constants : StefanBoltzmann_constant;
import nm.number;

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

import radiation.ray.tracing;

void main(string[] args) {

    // Default arg values
    string workingDir = ".";
    double absorptionCoefficient = 1.0;

    auto helpInformation = getopt(
        args, std.getopt.config.stopOnFirstNonOption,
        "d|dir", &workingDir,
        "k|absorptivity", &absorptionCoefficient
    );

    if (helpInformation.helpWanted) {
        defaultGetoptPrinter("lmr-raytrace options.", helpInformation.options);
        return;
    }

    writeln("Initial workings on a standalone radiation post-processing code.");
    writeln(format("Absorptivity: %.2g", absorptionCoefficient));

    // Hard-code how many snapshots we're working with
    uint snapshotStart = 0;
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

    writeln(format("Axisymmetric? %s", cfg.axisymmetric));

    // Set everything to zero initially
    // FIXME: This is because of some weird buffer thing in loading
    //        data during `lmr snapshot2vtk`, it takes the same
    //        values from the previous variable (temperature here)
    foreach (blk; localFluidBlocks) {
        foreach (cell; blk.cells) {
            cell.fs.Qrad = 0.0;
        }
    }

    string dirName = snapshotDirectory(nWrittenSnapshots);
    ensure_directory_is_present(dirName);

    FluidBlock block = localFluidBlocks[0];

    trace_intensity(block, absorptionCoefficient);

    foreach (blk; localFluidBlocks) {
        auto fileName = fluidFilename(nWrittenSnapshots, blk.id);
        FVCell[] cells;
        cells.length = blk.cells.length;
        foreach (i, ref c; cells) {
            c = blk.cells[i];
        }
        fluidBlkIO.writeVariablesToFile(fileName, cells);
    }
}
