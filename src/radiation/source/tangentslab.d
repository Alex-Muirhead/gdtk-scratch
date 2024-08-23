module tangentslab;

import std.math : exp, PI;

import nm.number;

/// A two term exponential approximation to the 3rd exponential
/// integral function E3(x)
///
/// Params:
///   x = Positive real input value
/// Returns: 
///   Approximate value
number exponential_integral(number x) {
    return 0.0929 * exp(-4.08 * x) + 0.4071 * exp(-1.33 * x);
}

number[] tangent_slab_heating(
    number[] intensity,
    number[] opticalThickness
) {
    size_t numCells = opticalThickness.length;
    number[] heating = new number[numCells];

    // Local variables
    number heatFlux;
    number halfThickness;
    number opticalDistance;

    foreach (i; 0 .. numCells) {
        // Do self-heating case here (i == j)
        halfThickness = opticalThickness[i] / 2;
        heatFlux = intensity[i] * (1 - 2 * exponential_integral(halfThickness));

        // Walk away from cell i (along negative coordinate)
        opticalDistance = halfThickness;
        for (auto j = long(i) - 1; j >= 0; --j) {
            heatFlux += intensity[j] * (
                exponential_integral(opticalDistance) -
                    exponential_integral(opticalDistance + opticalThickness[j])
            );
            opticalDistance += opticalThickness[j];
        }
        // Walk away from cell i (along positive coordinate)
        opticalDistance = halfThickness;
        for (auto j = i + 1; j < numCells; ++j) {
            heatFlux += intensity[j] * (
                exponential_integral(opticalDistance) -
                    exponential_integral(opticalDistance + opticalThickness[j])
            );
            opticalDistance += opticalThickness[j];
        }

        heating[i] = 2 * PI * heatFlux;
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
