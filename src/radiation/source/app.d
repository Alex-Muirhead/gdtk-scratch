import std.algorithm.mutation : swap;
import std.algorithm: min, max;
import std.array;
import std.math : cos, sin, PI, pow, exp, fabs;
import std.stdio;
import std.getopt;

import geom.elements.nomenclature : Face;
import gas.physical_constants : StefanBoltzmann_constant;
import geom.elements.vector3 : Vector3, distance_between;
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

import raytracing;


/// A two term exponential approximation to the 3rd exponential
/// integral function E3(x)
///
/// Params:
///   x = Positive real input value
/// Returns: 
///   Approximate value
number exponential_integral(number x) {
    return 0.0929*exp(-4.08*x) + 0.4071*exp(-1.33*x);
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

    string dirName = snapshotDirectory(nWrittenSnapshots);
    ensure_directory_is_present(dirName);

    double absorptionCoefficient = 0.5;
    FluidBlock block = localFluidBlocks[0];
    // trace_rays(block, decay);

    BoundaryCondition capsuleWall = block.bc[Face.east];
    size_t direction = Face.west;

    foreach (n, iface; capsuleWall.faces) {
        FluidFVCell cell = iface.left_cell;
        FluidFVCell[] tangentLine = [cell];

        // Keep track of the optical distances for the faces of each cell
        number opticalCoordinate = 0.0;
        number[] opticalCoordinates = [opticalCoordinate];

        number[] cellIntensity = [cell.grey_blackbody_intensity()];
        number[] cellHeatFlux = [];
        Vector3 lastPosition = cell.iface[Face.east].pos;
        Vector3 nextPosition;

        while (true) {
            nextPosition = cell.iface[direction].pos;
            // Absorption coefficient might change with temperature in the future
            opticalCoordinate += absorptionCoefficient * distance_between(nextPosition, lastPosition);
            opticalCoordinates ~= opticalCoordinate;

            if (cell.iface[direction].is_on_boundary) {
                // This appears to be the best way to check if we hit an edge
                // But it will also trigger if we cross a block boundary
                break;
            }
            // Cell cloud is offset by 1, as the cloud contains the cell itself
            cell = cell.cell_cloud[direction+1];

            tangentLine ~= cell;
            cellIntensity ~= cell.grey_blackbody_intensity();

            swap(nextPosition, lastPosition);
        }

        number heatFlux;
        number cellCenter;
        number edgeOne;
        number edgeTwo;
        number dropOff;
        
        foreach (i, ref c; tangentLine) {
            heatFlux = 0.0;
            cellCenter = (opticalCoordinates[i] + opticalCoordinates[i+1]) / 2;
            foreach (j, ref _c; tangentLine) {
                edgeOne = fabs(cellCenter - opticalCoordinates[j]);
                edgeTwo = fabs(cellCenter - opticalCoordinates[j+1]);
                if (j < i) {
                    // j+1 is closer than j
                    dropOff = exponential_integral(edgeTwo) - exponential_integral(edgeOne);
                } else if (j > i) {
                    // j is closer than j+1
                    dropOff = exponential_integral(edgeOne) - exponential_integral(edgeTwo);
                } else {
                    // In the same cell (i == j)       
                    dropOff = 0;
                }
                heatFlux += cellIntensity[j] * dropOff;
            }
            cellHeatFlux ~= heatFlux;
            c.fs.Qrad = heatFlux;

        }

        // DEBUG printing
        if (n < 5) {
            writeln("Optical distances from wall: ", opticalCoordinates);
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
