module raytracing;

// Standard modules
import std.conv : to;
import std.stdio;
import std.format;
import std.math;
import std.algorithm.mutation : swap;
import std.random : Random, uniform;

import gas.physical_constants : StefanBoltzmann_constant;
import geom.elements.vector3 : Vector3, wedge2D, dot;
import geom.grid.grid : Grid, Grid_t;
import nm.number : number;

import lmr.bc.boundary_condition : BoundaryCondition;
import lmr.bc.ghost_cell_effect.full_face_copy : GhostCellFullFaceCopy;
import lmr.fluidblock : FluidBlock;
import lmr.fluidfvcell : FluidFVCell;
import lmr.fvcell : FVCell;
import lmr.fvinterface : FVInterface;
import lmr.sfluidblock : SFluidBlock;
import lmr.ufluidblock : UFluidBlock;

import rays : Ray, HyperbolicRay, CartesianRay;

// External packages
import progress.bar;
import mir.random.engine;
import mirRandom = mir.random.engine;
import mir.random.ndvariable : sphereVar;
import mir.random.variable : uniformVar;

void trace_rays(FluidBlock block, number absorptionCoefficient) {
    auto rng = Random(4); // Chosen by fair dice roll guaranteed to be random (xkcd.com/221)
    // auto rne = mirRandom.Random(4);
    auto rne = mirRandom.Random(mirRandom.unpredictableSeed());
    uint angleSamples = 10_000;
    auto angleGenerator = uniformVar!number(0.0, 2 * PI);
    // NOTE: Using angleSamples of 1000 causes NaN values. Weird??

    number energyLost = 0.0;
    number energyEmitted = 0.0;
    number energyAbsorbed = 0.0;

    Bar progressBar = new Bar();
    progressBar.message = { return "Ray-tracing progress"; };
    progressBar.suffix = { return format("%6.2f%%\n", progressBar.percent); };
    progressBar.max = angleSamples * block.cells.length;

    number sign = 1.0;

    foreach (cell_id, ref origin; block.cells) {

        number fullEmissionEnergy = 4 * PI * origin.volume[0] * absorptionCoefficient * origin.grey_blackbody_intensity();
        writeln(format("[INFO ] Cell ID: %d,\t Energy: %.3g", cell_id, fullEmissionEnergy));
        // origin.fs.Qrad -= fullEmission; // Remove the energy of the ray emitted

        foreach (a; 0 .. angleSamples) {
            // number[3] angleVector;
            // sphereVar()(rne, angleVector); // TODO: Check if this is _true_ SO3

            // number angle = angleGenerator(rne);
            number angle = (double(a) / double(angleSamples) + 0.00001) * 2 * PI;

            // Vector3 direction = Vector3([angleVector[0], angleVector[1], 0]);
            // Vector3 direction = Vector3([angleVector[0], angleVector[1], angleVector[2]]);
            Vector3 direction = Vector3([0.0, cos(angle), sin(angle)]);

            // number mu = sqrt(1 - angleVector[2] ^^ 2);
            number mu = 1;

            number rayEnergy = fullEmissionEnergy / angleSamples;

            energyEmitted += rayEnergy;

            size_t[] rayCells;
            number[] rayLengths;

            FVInterface inter = marching_efficient(cell_id, block, direction, true, rayCells, rayLengths);

            //         writeln(format("[TRACE] Ray Cells: \t%s", rayCells));
            //         writeln(format("[TRACE] Ray lengths: \t%s", rayLengths));

            inner: foreach (i; 0 .. rayCells.length) {
                if (rayCells[i] == 1) {
                    cell_1_count += 1;
                    cell_1_length += rayLengths[i];
                }
                if (rayCells[i] == 4) {
                    cell_4_count += 1;
                    cell_4_length += rayLengths[i];
                }
                if (rayCells[i] == 7) {
                    cell_7_count += 1;
                    cell_7_length += rayLengths[i];
                }

                number cellVolume = block.cells[rayCells[i]].volume[0];
                // Kill off the ray if it's too weak
                if (rayEnergy < 1E-16) {
                    block.cells[rayCells[i]].fs.Qrad += rayEnergy / cellVolume;
                    break inner;
                }

                number opticalThickness = absorptionCoefficient * rayLengths[i] / mu;
                number heating = rayEnergy * (1 - exp(-opticalThickness));
                energyAbsorbed += heating;
                rayEnergy *= exp(-opticalThickness);
                block.cells[rayCells[i]].fs.Qrad += heating / cellVolume;
            }

            // writeln(format("Interface hit: %s", block.bc[inter.bc_id]));

            // Hopefully this reflects the ray properly?
            // FIXME: The new starting "position" isn't correct.
            // direction += -2 * direction.dot(inter.n) * inter.n;
            // direction.normalize();

            energyLost += rayEnergy;

            // progressBar.next();
        }
    }

    progressBar.finish();

    writeln(format("[DEBUG] Cell 1 count: %d, length: %.3f", cell_1_count, cell_1_length));
    writeln(format("[DEBUG] Cell 4 count: %d, length: %.3f", cell_4_count, cell_4_length));
    writeln(format("[DEBUG] Cell 7 count: %d, length: %.3f", cell_7_count, cell_7_length));

    writeln(format("Energy emitted: %.3g", energyEmitted));
    writeln(format("Energy absorbed: %.3g", energyAbsorbed));
    writeln(format("Energy lost: %.3g", energyLost));
    writeln(format("Discrepency: %.3g", energyEmitted - energyLost - energyAbsorbed));
}

FVInterface marching_efficient(
        size_t cellID, FluidBlock block, Vector3 rayTangent, bool axisymmetric,
        ref size_t[] crossed, ref number[] lengths
) {
    // Placeholder for if we need to do moving grids
    size_t gtl = 0;

    FluidFVCell currentCell = block.cells[cellID];
    Vector3 rayCoord = currentCell.pos[gtl];
    Grid currentGrid = get_grid(block);
    bool isWithinBlock = true;

    // Working variables
    FVInterface outgoing;
    Vector3 vertex_i;
    Vector3 vertex_j;

    rayTangent.normalize();

    // writeln(format("[INFO]  Tracing ray from %s in direction %s", rayCoord, rayTangent));

    Ray ray;
    if (axisymmetric) {
        ray = new HyperbolicRay(rayTangent, rayCoord);
    }
    else {
        ray = new CartesianRay(rayTangent, rayCoord);
    }

    // writeln(format("[DEBUG] Cell volume: %.3f", currentCell.volume[gtl]));
    // writeln(format("[INFO ] Starting at %s, in direction %s", rayCoord, rayTangent));
    // writeln(format("[INFO ] Ray description: %s", ray));

    while (isWithinBlock) {
    ndim:
        switch (currentGrid.dimensions) {
        case 1:
            throw new Exception("cell search not implemented for 1D grids");
        case 2:
            size_t iface_id;
            number step_length;

        faces: foreach (n, iface; currentCell.iface) {
                vertex_i = iface.vtx[0].pos[gtl];
                vertex_j = iface.vtx[1].pos[gtl];
                // Need to ensure the interface is moving anti-clockwise with 
                // respect to the cell interior
                if (currentCell.outsign[n] == -1) {
                    swap(vertex_i, vertex_j);
                }

                bool success = ray.intersect(vertex_i, vertex_j, step_length);
                if (success) {
                    iface_id = n;
                    break faces;
                }
            }

            debug {
                if (step_length.isNaN()) {
                    //                 writeln(format("[TRACE] Failed at cell %s", currentCell.id));
                    //                 writeln(format("[TRACE] Coords [%s, %s, %s, %s]", currentCell.vtx[0].pos, currentCell.vtx[1].pos, currentCell.vtx[2].pos, currentCell.vtx[3].pos));
                    auto hyper_ray = cast(HyperbolicRay) ray;
                    if (hyper_ray) {
                        //                     writeln(format("[TRACE] Failed with ray: \n\tdirection: %s", hyper_ray.localTangent(hyper_ray.linearCoord)));
                        //                     writeln(format("[TRACE] \tcoord: %s\n\ttangent: %s", rayCoord, rayTangent));
                    }

                    // Re-run this, and print debug values
                    foreach (n, iface; currentCell.iface) {
                        vertex_i = iface.vtx[0].pos[gtl];
                        vertex_j = iface.vtx[1].pos[gtl];
                        // Need to ensure the interface is moving anti-clockwise with 
                        // respect to the cell interior
                        if (currentCell.outsign[n] == -1) {
                            swap(vertex_i, vertex_j);
                        }

                        bool _ = ray.intersect(vertex_i, vertex_j, step_length);
                        //                     writeln(format("[TRACE] Checking vertex pair (%s, %s), got length %s", vertex_i, vertex_j, step_length));
                    }

                    throw new Exception("Something bad has happened...");
                    isWithinBlock = false;
                    break;
                }
            }

            ray.walkForward(step_length);
            outgoing = currentCell.iface[iface_id];

            lengths ~= step_length;
            crossed ~= currentCell.id;

            // debug {
            //     if ((crossed.length >= 4) && (crossed[$-1] == crossed[$-3]) && (crossed[$-2] == crossed[$-4])) {
            //         throw new Exception("Loop found within ray-tracing, exiting.");
            //     }
            // }

            // WARN: The `currentCell` can be null if crossing a boundary
            //       with no ghost cell on the other side
            currentCell = (currentCell.outsign[iface_id] == +1)
                ? outgoing.right_cell : outgoing.left_cell;

            // debug {
            //     writeln("Stepping ", step_length, " in direction ", rayTangent);
            //     writeln("Stepping across ", outgoing.id," to cell ", currentCell.id);
            // }

            if (outgoing.is_on_boundary) {
                isWithinBlock = false;
            }

            break ndim;
        case 3:
            throw new Exception("cell search not implemented for 3D grids");
        default:
            throw new Exception("invalid number of dimensions");
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
        size_t cellID, FluidBlock block, Vector3 rayTangent, bool axisymmetric,
        ref size_t[] crossed, ref number[] lengths
) {
    bool isWithinBlock = true;
    bool isContained;
    FluidFVCell currentCell = block.cells[cellID];
    Vector3 rayCoord = currentCell.pos[0]; // INFO: Grid is stationary, index is time
    Grid currentGrid = get_grid(block);

    Ray ray;
    if (axisymmetric) {
        ray = new HyperbolicRay(rayTangent, rayCoord);
    }
    else {
        ray = new CartesianRay(rayTangent, rayCoord);
    }

    number stepSize = 1E-03; // FIXME: This should be proportional to something...
    int consecutiveSteps = 1;
    FVInterface hit;

    // NOTE: Do we want to record the first step?

    marchLoop: while (isWithinBlock) {
        ray.walkForward(stepSize);
        rayCoord += stepSize * rayTangent;
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
                    hit = currentCell.iface[i - 1]; // 0th cell is self
                    break marchLoop;
                }
            }
            else {
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

void check_interface(FluidBlock block, FVInterface inter) {
    // Try and handle ghost cells
    BoundaryCondition boundary = block.bc[inter.bc_id];
    if (boundary.preReconAction.length > 0) {
        // writeln("Boundary type: ", boundary.type);
        auto mygce = cast(GhostCellFullFaceCopy) boundary.preReconAction[0];
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
            writeln("Interface ", inter.id, " maps to ", mappedCell.id);
        }
    }
}
