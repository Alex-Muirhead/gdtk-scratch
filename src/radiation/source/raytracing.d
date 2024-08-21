module raytracing;

import std.stdio;
import std.math;
import std.algorithm.mutation : swap;
import std.random : Random, uniform;

import gas.physical_constants : StefanBoltzmann_constant;
import geom.elements.vector3 : Vector3, wedge2D;
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

void trace_rays(FluidBlock block, number absorptivity) {
    auto rng = Random(4); // Chosen by fair dice roll guaranteed to be random (xkcd.com/221)
    uint angleSamples = 100;

    // FIXME: This needs to scan over for the minimum temperature
    number baseEmission = StefanBoltzmann_constant * 280 ^^ 4;

    outer: foreach (cell_id, ref origin; block.cells) {

        // auto angle = uniform(0.0f, 2*PI, rng);
        number fullEmission = origin.volume[0] * StefanBoltzmann_constant * origin.fs.gas.T ^^ 4;
        fullEmission -= baseEmission; // Add back baseline background radiation
        origin.fs.Qrad -= fullEmission; // Remove the energy of the ray emitted

        foreach (a; 0 .. angleSamples) {
            // Loop through angles for now
            auto angle = 2 * PI * a / angleSamples;
            auto direction = Vector3([cos(angle), sin(angle), 0]);
            direction.normalize(); // Shouldn't be needed with random angle

            auto rayStrength = fullEmission / angleSamples;

            size_t[] rayCells;
            number[] rayLengths;
            // FVInterface inter = marching_full(cell_id, firstBlock, direction, rayCells, rayLengths);
            FVInterface inter = marching_efficient(cell_id, block, direction, rayCells, rayLengths);

            // Try and handle ghost cells
            BoundaryCondition boundary = block.bc[inter.bc_id];
            if (boundary.preReconAction.length > 0) {
                writeln("Boundary type: ", boundary.type);
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
                    writeln("Cell ", rayCells[$ - 1], " maps to ", mappedCell.id);
                }
            }

            number heating;
            foreach (i; 0 .. rayCells.length) {
                heating = rayStrength * (1 - (1 - absorptivity) ^^ rayLengths[i]);
                rayStrength -= heating;
                block.cells[rayCells[i]].fs.Qrad += heating;
            }
        }
        break outer;
    }
}

FVInterface marching_efficient(
    size_t cellID, FluidBlock block, Vector3 rayTangent,
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

    // NOTE: Push all this down into the grid level
    //       Potentially with the following signature

    /**
     * Params:
     *   line = A description of the ray position & orientation
     *   indx = The cell index to check for intersections at (note, is this enough info?)
     *   crossings = The interfaces which the line crosses through. The two (2) solutions correspond
     *               to the nonpositive (<=0) and positive (>0) solutions.
     */
    // bool compute_line_intersections(ref Ray line, ref size_t indx, out FVInterface[2] crossings)

    rayTangent.normalize();

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

                Vector3 toVertex = rayCoord - vertex_i;
                // NOTE: Could potentially use iface.t1 for this?
                //       Probably not, since we need a non-normalised vec 
                Vector3 faceTangent = vertex_j - vertex_i;

                number alignment = wedge2D(faceTangent, rayTangent);
                number s = wedge2D(toVertex, faceTangent) / alignment;
                number t = wedge2D(toVertex, rayTangent) / alignment;

                // There should only be a single solution if the cell is convex
                if (t >= 0 && t < 1 && s > number.epsilon) {
                    iface_id = n;
                    step_length = s;
                    break faces;
                }
            }

            rayCoord += rayTangent * step_length;
            outgoing = currentCell.iface[iface_id];

            lengths ~= step_length;
            crossed ~= currentCell.id;

            // WARN: The `currentCell` can be null if crossing a boundary
            //       with no ghost cell on the other side
            currentCell = (currentCell.outsign[iface_id] == +1)
                ? outgoing.right_cell : outgoing.left_cell;

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
