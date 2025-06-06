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
import lmr.fvcell : FVCell;
import lmr.fvinterface : FVInterface;
import lmr.globalconfig;
import lmr.globaldata;
import lmr.init;
import lmr.sfluidblock;
import lmr.ufluidblock;

import singlog : logger = log;

import radiation.ray.tracing;

void main(string[] args) {

    // Default arg values
    string workingDir = ".";
    double absorptionCoefficient = 1.0;
    bool modelEmission = false;

    auto helpInformation = getopt(
        args, std.getopt.config.stopOnFirstNonOption,
        "d|dir", &workingDir,
        "k|absorptivity", &absorptionCoefficient,
        "emission", &modelEmission
    );

    if (helpInformation.helpWanted) {
        defaultGetoptPrinter("lmr-raytrace options.", helpInformation.options);
        return;
    }

    logger.program("Ray Tracing")
        .color(true)
        // .level(logger.DEBUGGING)
        .output(logger.output.std);

    writeln("Initial workings on a standalone radiation post-processing code.");
    logger.information(format("Absorptivity: %.2g", absorptionCoefficient));

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

    logger.information(format("Axisymmetric? %s", cfg.axisymmetric));

    // Set initial value before tracing
    // Can use the default emission for simplicity
    foreach (blk; localFluidBlocks) {
        foreach (cell; blk.cells) {
            if (modelEmission) {
                cell.fs.Qrad = -4*PI* absorptionCoefficient * cell.grey_blackbody_intensity();
            } else {
                cell.fs.Qrad = 0.0;
            }
        }
    }

    string dirName = snapshotDirectory(nWrittenSnapshots);
    ensure_directory_is_present(dirName);

    FluidBlock block = localFluidBlocks[0];
    trace_rays(block, absorptionCoefficient);

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
