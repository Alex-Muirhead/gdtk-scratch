module radiation.ray.tracing;

// Standard modules
import std.algorithm.iteration : sum, map;
import std.conv : to;
import std.stdio;
import std.format;
import std.math;
import std.algorithm.mutation : swap;
import std.random : Random, uniform;

import gas.physical_constants : StefanBoltzmann_constant;
import geom.elements.vector3 : Vector3, wedge2D, dot;
import geom.elements.properties : on_left_of_xy_line;
import geom.grid.grid : Grid, Grid_t;
import nm.number : number;

import lmr.bc.boundary_condition : BoundaryCondition;
import lmr.bc.ghost_cell_effect.full_face_copy : GhostCellFullFaceCopy;
import lmr.fluidblock : FluidBlock;
import lmr.fluidfvcell : FluidFVCell;
import lmr.fvcell : FVCell;
import lmr.fvinterface : FVInterface;
import lmr.globalconfig : GlobalConfig;
import lmr.sfluidblock : SFluidBlock;
import lmr.ufluidblock : UFluidBlock;

alias cfg = GlobalConfig;

import radiation.ray.types : Ray, HyperbolicRay, PlanarRay;

// External packages
import mir.random.engine;
import mirRandom = mir.random.engine;
import mir.random.ndvariable : sphereVar;
import mir.random.variable : uniformVar;
import progress.bar;
import singlog : logger = log;  // `log` is already bound to `math.log`

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

    foreach (cell_id, ref origin; block.cells) {

        number fullEmissionEnergy = 4 * PI * origin.volume[0]
            * absorptionCoefficient * origin.grey_blackbody_intensity();
        // writeln(format("[INFO ] Cell ID: %d,\t Energy: %.3g", cell_id, fullEmissionEnergy));

        number selfAbsorbed = 0.0;

        foreach (a; 0 .. angleSamples) {
            number[3] angleVector;
            sphereVar()(rne, angleVector); // TODO: Check if this is _true_ SO3

            // number angle = angleGenerator(rne);
            // number angle = (double(a) / double(angleSamples) + 1E-6) * 2 * PI;

            // Vector3 direction = Vector3([angleVector[0], angleVector[1], 0]);
            Vector3 direction = Vector3([angleVector[0], angleVector[1], angleVector[2]]);
            // Vector3 direction = Vector3([0.0, cos(angle), sin(angle)]);

            // number mu = sqrt(1 - angleVector[2] ^^ 2);
            number mu = 1;

            number rayEnergy = fullEmissionEnergy / angleSamples;

            energyEmitted += rayEnergy;

            size_t[] rayCells;
            number[] rayLengths;

            FVInterface inter = marching_efficient(cell_id, block, direction, cfg.axisymmetric, rayCells, rayLengths);

            number initialEnergy = rayEnergy;
            inner: foreach (i; 0 .. rayCells.length) {
                number cellVolume = block.cells[rayCells[i]].volume[0];
                // Kill off the ray if it's too weak
                if (rayEnergy / initialEnergy < 1E-16) {
                    block.cells[rayCells[i]].fs.Qrad += rayEnergy / cellVolume;
                    energyAbsorbed += rayEnergy;
                    rayEnergy = 0.0;
                    break inner;
                }

                number opticalThickness = absorptionCoefficient * rayLengths[i] / mu;
                number heating = rayEnergy * (1 - exp(-opticalThickness));
                if (rayCells[i] == cell_id) {
                    selfAbsorbed += heating;
                }
                energyAbsorbed += heating;
                rayEnergy *= exp(-opticalThickness);
                block.cells[rayCells[i]].fs.Qrad += heating / cellVolume;
            }

            // writeln(format("Interface hit: %s", block.bc[inter.bc_id]));

            energyLost += rayEnergy;

            progressBar.next();
        } // foreach (a; 0 .. angleSamples)
    }

    progressBar.finish();

    writeln(format("Energy emitted: %.3g", energyEmitted));
    writeln(format("Energy absorbed: %.3g", energyAbsorbed));
    writeln(format("Energy lost: %.3g", energyLost));
    writeln(format("Discrepency: %.3g", energyEmitted - energyLost - energyAbsorbed));
}

void trace_intensity(FluidBlock block, number absorptionCoefficient) {
    // auto rng = Random(4); // Chosen by fair dice roll guaranteed to be random (xkcd.com/221)
    // auto rne = mirRandom.Random(4);
    // auto rne = mirRandom.Random(mirRandom.unpredictableSeed());
    uint angleSamples = 1_000;
    number anglePatch = 4 * PI / double(angleSamples);
    // auto angleGenerator = uniformVar!number(0.0, 2 * PI);

    Bar progressBar = new Bar();
    progressBar.message = { return "Ray-tracing progress"; };
    progressBar.suffix = { return format("%6.2f%%\n", progressBar.percent); };
    progressBar.max = angleSamples * block.cells.length;

    number expectedIntensity = block.cells[0].grey_blackbody_intensity();
    writeln(format("Expected intensity: %.3g", expectedIntensity));

    foreach (cell_id, ref origin; block.cells) {

        number incidentRadiation = 0.0;

        foreach (a; 0 .. angleSamples) {
            number angle = (double(a) / double(angleSamples) + 1E-6) * 2 * PI;
            Vector3 direction = Vector3([0.0, cos(angle), sin(angle)]);

            size_t[] rayCells;
            number[] rayLengths;

            FVInterface _ = marching_efficient(cell_id, block, direction, true, rayCells, rayLengths);
            // FVInterface inter = marching_full(cell_id, block, direction, true, rayCells, rayLengths);

            // This is entirely contained within the ray
            number opticalDistance = 0.0;
            foreach (i; 0 .. rayCells.length) {
                number opticalThickness = absorptionCoefficient * rayLengths[i];
                number uniformIntensity = block.cells[rayCells[i]].grey_blackbody_intensity();
                number exitIntensity = uniformIntensity * (1 - exp(-opticalThickness));
                number impactedIntensity = exitIntensity * exp(-opticalDistance);
                // G = integral vec(I) dot dd(vec(Omega))
                // vec(q) = integral vec(I) dd(Omega)
                incidentRadiation += impactedIntensity * anglePatch;
                opticalDistance += opticalThickness;
            }

            progressBar.next();
        }

        origin.fs.Qrad += absorptionCoefficient * incidentRadiation;
    }

    progressBar.finish();
}

struct Crossing {
    size_t iface_id;
    number distance;
}

void trace_block(FluidBlock block, FluidFVCell starting_cell, Ray ray) {

}

bool trace_cell(FluidFVCell cell, Ray ray, out size_t iface_id, out number step_length) {
    // Placeholder for if we need to do moving grids
    size_t gtl = 0;

    Vector3 vertex_i;
    Vector3 vertex_j;

    debug {
        logger.debugging(format("Coords %s", cell.vtx.map!(v => v.pos[gtl])));
    }

    foreach (n, iface; cell.iface) {
        vertex_i = iface.vtx[0].pos[gtl];
        vertex_j = iface.vtx[1].pos[gtl];
        // Need to ensure the interface is moving anti-clockwise with
        // respect to the cell interior
        if (cell.outsign[n] == -1) {
            swap(vertex_i, vertex_j);
        }

        bool success = ray.intersect(vertex_i, vertex_j, step_length);
        debug {
            logger.debugging(format(
                    "Checking vertex pair (%s, %s), got length %s",
                    vertex_i, vertex_j, step_length
            ));
        }
        if (success) {
            iface_id = n;
            return true;
        }
    }

    return false;
}

FVInterface marching_efficient(
    size_t cellID, FluidBlock block,
    Vector3 rayTangent, bool axisymmetric,
    ref size_t[] crossed, ref number[] lengths
) {
    // Placeholder for if we need to do moving grids
    size_t gtl = 0;

    FluidFVCell currentCell = block.cells[cellID];
    Vector3 rayCoord = currentCell.pos[gtl];

    rayTangent.normalize();

    // TODO: Should we push this to a generic, to specialise?
    Ray ray;
    if (axisymmetric) {
        ray = new HyperbolicRay(rayTangent, rayCoord);
    } else {
        ray = new PlanarRay(rayTangent, rayCoord);
    }

    FVInterface outgoing;

    bool isWithinBlock = true;
    while (isWithinBlock) {
        size_t iface_id;
        number step_length;

        bool success = trace_cell(currentCell, ray, iface_id, step_length);
        debug {
            // Re-run if necessary
            if (!success) {
                logger.level(logger.DEBUGGING);
                logger.error(format("Failed at cell %s. Re-running with verbosity.", currentCell.id));
                logger.debugging(format("Ray: %s", ray));
                auto _ = trace_cell(currentCell, ray, iface_id, step_length);
                logger.level(logger.ERROR);
            }
        }

        // Does this actually help us much? Or should we just leave the ray as static,
        // and solve both intersections for each cell?
        ray.walkForward(step_length);
        outgoing = currentCell.iface[iface_id];

        lengths ~= step_length;
        crossed ~= currentCell.id;

        // WARN: The `currentCell` can be null if crossing a boundary
        //       with no ghost cell on the other side
        currentCell = (currentCell.outsign[iface_id] == +1) ? outgoing.right_cell : outgoing.left_cell;

        // What other terminating conditions might we want?
        if (outgoing.is_on_boundary) {
            isWithinBlock = false;
        }

    }

    return outgoing;
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
