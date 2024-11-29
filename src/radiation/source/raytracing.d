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

// External packages
import progress.bar;
import mir.random.engine;
import mirRandom = mir.random.engine;
import mir.random.ndvariable : sphereVar;

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

void trace_rays(FluidBlock block, number absorptionCoefficient) {
    auto rng = Random(4); // Chosen by fair dice roll guaranteed to be random (xkcd.com/221)
    auto rne = mirRandom.Random(4);
    uint angleSamples = 967;
    // NOTE: Using angleSamples of 1000 causes NaN values. Weird??

    number energyLost = 0.0;
    number energyEmitted = 0.0;
    number energyAbsorbed = 0.0;

    Bar progressBar = new Bar();
    progressBar.message = { return "Ray-tracing progress"; };
    progressBar.suffix = { return format("%6.2f%%", progressBar.percent); };
    progressBar.max = angleSamples * block.cells.length;

    foreach (cell_id, ref origin; block.cells) {

        number fullEmissionEnergy = 4 * PI * origin.volume[0] * absorptionCoefficient * origin.grey_blackbody_intensity();
        // origin.fs.Qrad -= fullEmission; // Remove the energy of the ray emitted

        foreach (a; 0 .. angleSamples) {
            number[3] angleVector;
            sphereVar()(rne, angleVector); // TODO: Check if this is _true_ SO3

            // Vector3 direction = Vector3([angleVector[0], angleVector[1], 0]);
            Vector3 direction = Vector3([angleVector[0], angleVector[1], angleVector[2]]);
            // number mu = sqrt(1 - angleVector[2] ^^ 2);
            number mu = 1;

            number rayEnergy = fullEmissionEnergy / angleSamples;

            energyEmitted += rayEnergy;

            size_t[] rayCells;
            number[] rayLengths;
            FVInterface inter = marching_efficient(cell_id, block, direction, true, rayCells, rayLengths);

            inner: foreach (i; 0 .. rayCells.length) {
                number cellVolume = block.cells[rayCells[i]].volume[0];
                // Kill off the ray if it's too weak
                if (rayEnergy < 1E-10) {
                    block.cells[rayCells[i]].fs.Qrad += rayEnergy / cellVolume;
                    break inner;
                }

                number opticalThickness = absorptionCoefficient * rayLengths[i] / mu;
                number heating = rayEnergy * (1 - exp(-opticalThickness));
                energyAbsorbed += heating;
                rayEnergy *= exp(-opticalThickness);
                block.cells[rayCells[i]].fs.Qrad += heating / cellVolume;
            }
            energyLost += rayEnergy;

            progressBar.next();
        }
    }

    progressBar.finish();

    writeln(format("Energy emitted: %.3g", energyEmitted));
    writeln(format("Energy absorbed: %.3g", energyAbsorbed));
    writeln(format("Energy lost: %.3g", energyLost));
    writeln(format("Discrepency: %.3g", energyEmitted - energyLost - energyAbsorbed));
}

interface Ray {
    void walkForward(number length);
    bool intersect(Vector3 vertexOne, Vector3 vertexTwo, out number length);
}

class CartesianRay : Ray {
public:
    Vector3 terminus;
    Vector3 tangent;
    number coord;

    this(Vector3 tangent, Vector3 terminus) {
        this.tangent = tangent;
        this.terminus = terminus;
    }

    /**
     * Moves the ray forward by the specified length.
     *
     * Params:
     *   length = The distance to move forward (number).
     */
    void walkForward(number length) {
        coord += length;
    }

    /**
     * Checks if the ray intersects the line segment between two vertices.
     *
     * Params:
     *   vertexOne = The first vertex of the line segment (Vector3).
     *   vertexTwo = The second vertex of the line segment (Vector3).
     *   length = The arc length at the intersection point (out parameter) (number).
     *
     * Returns:
     *   true if the ray intersects the line segment, false otherwise.
     */
    bool intersect(Vector3 vertexOne, Vector3 vertexTwo, out number length) {
        // Calculate the face tangent vector
        // NOTE: Could potentially use iface.t1 for this?
        //       Probably not, since we need a non-normalised vec 
        Vector3 faceTangent = vertexTwo - vertexOne;
        Vector3 toVertex = terminus - vertexOne;

        number alignment = wedge2D(faceTangent, tangent);
        number s = wedge2D(toVertex, faceTangent) / alignment;
        number t = wedge2D(toVertex, tangent) / alignment;

        // There should only be a single solution if the cell is convex
        if (t >= 0 && t < 1 && s > 1E5 * number.epsilon) {
            length = s;
            return true;
        }

        return false;
    }
}

class HyperbolicRay : Ray {
public:
    Vector3 center;
    number radius;
    number slope;
    number arcCoord;
private:
    Vector3 positiveAsymptote;
    Vector3 negativeAsymptote;
    number correction;

public:
    /**
     * Constructs a HyperbolicRay instance.
     *
     * Params:
     *   tangent = The tangent vector (Vector3).
     *   point = The point on the ray (Vector3).
     */
    this(Vector3 tangent, Vector3 terminus) {
        arcCoord = (terminus.y * tangent.y + terminus.z * tangent.z)
            / (tangent.y ^^ 2 + tangent.z ^^ 2);

        Vector3 saddle = terminus - tangent * arcCoord;

        center = Vector3(saddle.x, 0.0);
        radius = sqrt(saddle.y ^^ 2 + saddle.z ^^ 2);
        slope = tangent.x / sqrt(1 - tangent.x ^^ 2);

        // Positive and Negative are defined relative to tangent vector.
        // The arc-coordinate is defined in the same frame.
        positiveAsymptote = Vector3(+slope, 1);
        negativeAsymptote = Vector3(-slope, 1);

        // writeln(format("[DEBUG] Hyperbolic ray -> Center: %s, Radius: %s, Slope: %f", center, radius, slope));
        // writeln(format("                       -> Asymptote pos: %s, neg: %s", positiveAsymptote, negativeAsymptote));

        correction = sqrt(1 + slope ^^ 2); // Magnitude of asymptote vectors
    }

    Vector3 localTangent(number localArcCoord) {
        // (positiveAsymptote - negativeAsymptote) / 2
        // and (positiveAsymptote + negativeAsymptote) / 2
        return Vector3(slope, 0.0) + Vector3(0.0, 1.0) * localArcCoord / correction 
            / sqrt((localArcCoord / correction) ^^ 2 + radius ^^ 2);
    }

    /**
     * Moves the ray forward by the specified arc length.
     *
     * Params:
     *   arcLength = The distance to move forward (number).
     */
    void walkForward(number arcLength) {
        arcCoord += arcLength;
    }

    /**
     * Checks if the ray intersects the line segment between two vertices.
     *
     * Params:
     *   vertexOne = The first vertex of the line segment (Vector3).
     *   vertexTwo = The second vertex of the line segment (Vector3).
     *   arcLength = The arc length at the intersection point (out parameter) (number).
     *
     * Returns:
     *   true if the ray intersects the line segment, false otherwise.
     */
    bool intersect(Vector3 vertexOne, Vector3 vertexTwo, out number arcLength) {
        // Calculate the face tangent vector
        Vector3 faceTangent = vertexTwo - vertexOne;

        // Calculate the alignment of the face tangent with the asymptotes
        number positiveAlignment = wedge2D(positiveAsymptote, faceTangent);
        number negativeAlignment = wedge2D(negativeAsymptote, faceTangent);

        // writeln(format("[DEBUG] Positive Alignment: %e", positiveAlignment));
        // writeln(format("[DEBUG] Negative Alignment: %e", negativeAlignment));

        // Calculate the vector from the center to the first vertex
        Vector3 toVertex = vertexOne - center;
        number findAName = wedge2D(toVertex, faceTangent);
        number determinant = findAName ^^ 2
            - radius ^^ 2 * positiveAlignment * negativeAlignment;

        // writeln(format("[DEBUG] Determinant: %e", determinant));
        // Check if the determinant is negative (no intersection)
        if (determinant < 0) {
    //         writeln("[INFO]  Skipped face");
            return false;
        }

        // Try both solutions using static foreach
        // NOTE: `static foreach` has the limitation that you can't use `continue` etc
        // OH WAIT! We only use the `sign` in the two exponential solutions, but we can
        // Throw away one of them entirely
        // ...but there might be two solutions still if sqrt(determinant) < findAName
        foreach (sign; [+1, -1]) {
            {
                // Calculate the exponential, arc, Cartesian, and linear coordinates at the intercept
                number exponential = (findAName + sign * sqrt(determinant)) / positiveAlignment;
                if (exponential < 0) {
                //     writeln(format("[TRACE] Skipping mirror solution"));
                    continue;
                }

                number arc = (exponential ^^ 2 - radius ^^ 2) / (2 * exponential) * correction;
                number crossingDirection = wedge2D(localTangent(arc), faceTangent);
        //         writeln(format("[DEBUG] Local tangent: %s", localTangent(arc)));
        //         writeln(format("[DEBUG] Crossing direction: %e", crossingDirection));

                Vector3 cartesian = center
                    + positiveAsymptote * (exponential / 2)
                    + negativeAsymptote / (exponential * 2) * radius ^^ 2;
                Vector3 translation = cartesian - vertexOne;
                number linear = dot(translation, faceTangent) / dot(faceTangent, faceTangent);

        //         writeln(format("[DEBUG] Pos part: %s", positiveAsymptote * (exponential / 2)));
        //         writeln(format("[DEBUG] Neg part: %s", negativeAsymptote / (exponential * 2) * radius ^^ 2));

        //         writeln(format("[DEBUG] Exponential coord: %f", exponential));
        //         writeln(format("[DEBUG] Arc step: %e, from: %e, to: %s", arc - arcCoord, arcCoord, arc));
        //         writeln(format("[DEBUG] Cartesian coord: %s", cartesian));
        //         writeln(format("[DEBUG] Linear coord: %f", linear));

                // Check if the intersection point is within the line segment
                // and ahead of the current ray position
                if (
                    linear >= 0 && linear < 1
                    // FIXME: Is there a more robust way of checking?
                    // i.e. comparing the current tangent to the face tangent?
                    // since we only care about *leaving* the current face.
                    // We could add a parameter for positive or negative crossing
                    && arc > arcCoord
                    && crossingDirection > 0
                ) {
                    arcLength = arc - arcCoord;
                    return true;
                }
            }
        }

        return false;
    }
}

@("Hyperbolic Intercept - Triangle")
unittest {
    import fluent.asserts;

    Vector3 point = Vector3(2.0, 1.5, 0.0);
    Vector3 tangent = Vector3(-1.5, -0.5, 1.5);
    tangent.normalize();

    // writeln(format("Starting direction: %s", tangent));

    Ray ray = new HyperbolicRay(tangent, point);

    // Define triangle vertices
    Vector3 vertex0 = Vector3(0.25, 0.5);
    Vector3 vertex1 = Vector3(3.25, 1.25);
    Vector3 vertex2 = Vector3(1.5, 3.0);

    number arcLength;

    // Test intersection with triangle edge vertex2 to vertex0
    bool success = ray.intersect(vertex2, vertex0, arcLength);
    // writeln(format("Arclength: %.7f", arcLength));
    Assert.equal(success, true);
    Assert.approximately(arcLength, 1.72741, 1e-5);

    // Test intersection with triangle edge vertex1 to vertex2 (should fail)
    success = ray.intersect(vertex1, vertex2, arcLength);
    Assert.equal(success, false);
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
    } else {
        ray = new CartesianRay(rayTangent, rayCoord);
    }

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

                // writeln(format("[TRACE] Checking vertex pair %s, %s", vertex_i, vertex_j));

                bool success = ray.intersect(vertex_i, vertex_j, step_length);
                if (success) {
    //                 writeln(format("[INFO]  Crossing interface %d", n));
                    iface_id = n;
                    break faces;
                }
            }

            if (step_length.isNaN()) {
                throw new Exception("Something bad has happened...");
                isWithinBlock = false;
                break;
            }

            ray.walkForward(step_length);
            outgoing = currentCell.iface[iface_id];

            lengths ~= step_length;
            crossed ~= currentCell.id;

            // WARN: The `currentCell` can be null if crossing a boundary
            //       with no ghost cell on the other side
            currentCell = (currentCell.outsign[iface_id] == +1)
                ? outgoing.right_cell : outgoing.left_cell;
            // writeln("Stepping ", step_length, " in direction ", rayTangent);
            // writeln("Stepping across ", outgoing," to cell ", currentCell);

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
    size_t cellID, FluidBlock block, Vector3 direction,
    ref size_t[] crossed, ref number[] lengths) {
    bool isWithinBlock = true;
    bool isContained;
    FluidFVCell currentCell = block.cells[cellID];
    Vector3 rayCoord = currentCell.pos[0]; // INFO: Grid is stationary, index is time
    Grid currentGrid = get_grid(block);

    number stepSize = 1E-03; // FIXME: This should be proportional to something...
    int consecutiveSteps = 1;
    FVInterface hit;

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
