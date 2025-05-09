import std.algorithm.mutation : swap;
import std.algorithm : min, max;
import std.array;
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

import tangentslab;

void main(string[] args) {

    string workingDir = ".";
    getopt(args, std.getopt.config.stopOnFirstNonOption, "d|dir", &workingDir);

    writeln("Initial workings on a standalone radiation post-processing code.");

    // Hard-code how many snapshots we're working with
    uint snapshotStart = 0000;
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

    string dirName = snapshotDirectory(nWrittenSnapshots);
    ensure_directory_is_present(dirName);

    double absorptionCoefficient = 4.0;
    FluidBlock block = localFluidBlocks[0];

    uint wallSide = Face.east;
    uint direction = opposite_face(wallSide);
    BoundaryCondition capsuleWall = block.bc[wallSide];

    foreach (n, iface; capsuleWall.faces) {

        // Start tangent trace routine

        size_t[] crossed = [];
        number[] lengths = [];

        naive_tangent_marching(iface.left_cell.id, block, direction, crossed, lengths);
        // NOTE: Need to switch interface SIDE and normal SIGN when going from West to East etc.
        // marching_efficient(iface.right_cell.id, block, iface.n, crossed, lengths);

        number[] cellIntensity = [];
        number[] opticalThickness = [];

        for (auto i = 0; i < crossed.length; i++) {
            // NOTE: Absorption coefficient might be temperature dependent
            opticalThickness ~= absorptionCoefficient * lengths[i];
            cellIntensity ~= block.cells[crossed[i]].grey_blackbody_intensity();
        }

        number[] cellHeating = tangent_slab_heating(cellIntensity, opticalThickness);

        foreach (i, cellID; crossed) {
            block.cells[cellID].fs.Qrad = cellHeating[i] * block.cells[cellID].volume[0];
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

FVInterface naive_tangent_marching(size_t cellID, FluidBlock block, uint rayDirection,
    ref size_t[] crossed, ref number[] lengths) {
    FluidFVCell currentCell = block.cells[cellID];
    // Grid currentGrid = get_grid(block); // NOTE: Could use the grid for direction?

    while (true) {
        crossed ~= currentCell.id;
        lengths ~= distance_between(currentCell.iface[rayDirection].pos,
            currentCell.iface[opposite_face(rayDirection)].pos);

        // This appears to be the best way to check if we hit an edge
        // But it will also trigger if we cross a block boundary
        if (currentCell.iface[rayDirection].is_on_boundary) {
            break;
        }

        // Cell cloud is offset by 1, as the cloud contains the cell itself
        currentCell = currentCell.cell_cloud[rayDirection + 1];
    }

    return currentCell.iface[rayDirection];
}
