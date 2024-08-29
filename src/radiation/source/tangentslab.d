module tangentslab;

import std.math : exp, PI;
import std.stdio;
import std.format;

import nm.number;

/// A two term exponential approximation to the 3rd exponential
/// integral function E_2(x)
///
/// Params:
///   x = Positive real input value
/// Returns: 
///   Approximate value
number exponential_integral_2(number x) {
    return 0.6930 * exp(-1.52 * x) + 0.3069 * exp(-8.45 * x);
}

/// A two term exponential approximation to the 3rd exponential
/// integral function E_3(x)
///
/// Params:
///   x = Positive real input value
/// Returns: 
///   Approximate value
number exponential_integral_3(number x) {
    return 0.0929 * exp(-4.08 * x) + 0.4071 * exp(-1.33 * x);
}

number[] tangent_slab_flux(
    number[] intensity,
    number[] opticalThickness
) {
    size_t numCells = opticalThickness.length;
    number[] heating = new number[numCells];

    // Local variables
    number heatFlux;
    number halfThickness;
    number opticalDistance;

    number increment;

    foreach (i; 0 .. numCells) {
        // Do self-heating case here (i == j)
        halfThickness = opticalThickness[i] / 2;
        heatFlux = 0;

        // Walk away from cell i (along negative coordinate)
        opticalDistance = halfThickness;
        for (auto j = long(i) - 1; j >= 0; --j) {
            increment = intensity[j] * (
                exponential_integral_3(opticalDistance) -
                    exponential_integral_3(opticalDistance + opticalThickness[j])
            );
            heatFlux += increment;
            writeln(format("Cell %d increases %d by %f", j, i, increment));
            opticalDistance += opticalThickness[j];
        }
        // Walk away from cell i (along positive coordinate)
        opticalDistance = halfThickness;
        for (auto j = i + 1; j < numCells; ++j) {
            increment = intensity[j] * (
                exponential_integral_3(opticalDistance + opticalThickness[j]) - 
                    exponential_integral_3(opticalDistance)
            );
            heatFlux += increment;
            writeln(format("Cell %d increases %d by %f", j, i, increment));
            opticalDistance += opticalThickness[j];
        }

        heating[i] = 2 * PI * heatFlux;
    }

    return heating;
}

number[] tangent_slab_heating(
    number[] intensity,
    number[] opticalThickness
) {
    size_t numCells = opticalThickness.length;
    number[] heating = new number[numCells];

    // Local variables
    number incidentRadiation;
    number halfThickness;
    number opticalDistance;

    number increment;

    foreach (i; 0 .. numCells) {
        writeln(format("Cell %d intensity: %f", i, intensity[i]));
        // Do self-heating case here (i == j)
        halfThickness = opticalThickness[i] / 2;
        incidentRadiation = intensity[i] * (2 - 2 * exponential_integral_2(halfThickness));
        writeln(format("Cell %d increases %d by %f", i, i, incidentRadiation));

        // Walk away from cell i (along negative coordinate)
        opticalDistance = halfThickness;
        for (auto j = long(i) - 1; j >= 0; --j) {
            increment = intensity[j] * (
                exponential_integral_2(opticalDistance) -
                    exponential_integral_2(opticalDistance + opticalThickness[j])
            );
            incidentRadiation += increment;
            writeln(format("Cell %d increases %d by %f", j, i, increment));
            opticalDistance += opticalThickness[j];
        }
        // Walk away from cell i (along positive coordinate)
        opticalDistance = halfThickness;
        for (auto j = i + 1; j < numCells; ++j) {
            increment = intensity[j] * (
                exponential_integral_2(opticalDistance) -
                    exponential_integral_2(opticalDistance + opticalThickness[j])
            );
            incidentRadiation += increment;
            writeln(format("Cell %d increases %d by %f", j, i, increment));
            opticalDistance += opticalThickness[j];
        }

        heating[i] = 2 * PI * incidentRadiation;
    }

    return heating;
}

@("Constant Temperature")
unittest {

    import std.stdio;
    import std.range;
    import std.algorithm.iteration;
    import std.math.operations;

    const size_t numCells = 10;

    number[numCells] intensity = 1.0;
    number[numCells] opticalThickness = 1.0 / numCells;
    number[numCells] cellCenters;

    cellCenters[0] = opticalThickness[0] / 2;
    for (auto i = 1; i < numCells; ++i) {
        cellCenters[i] = cellCenters[i-1] + (opticalThickness[i] + opticalThickness[i-1]) / 2;
    }

    // D is so weird, wtf do I need to slice a static array to sum it
    number opticalLength = opticalThickness[].sum();

    number[numCells] heating = tangent_slab_heating(intensity, opticalThickness);

    for (auto i = 0; i < numCells; ++i) {
        // Analytical solution
        number expected = 2 * PI * (
            2 * exponential_integral(0) - 
            exponential_integral(cellCenters[i]) - 
            exponential_integral(opticalLength - cellCenters[i])
        );
        assert(isClose(heating[i], expected, 1e-6, 1e-6));
    }
}
