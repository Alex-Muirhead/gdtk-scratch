module tangentslab;

import std.math : exp;

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
    number[] cellIntensity,
    number[] opticalThickness,
    number[] lengths
) {
    size_t numCells = lengths.length;
    number[] heating;
    heating.length = numCells;

    // Local variables
    number heatFlux;
    number halfThickness;
    number opticalDistance;

    foreach (i; 0 .. numCells) {
        // Do self-heating case here (i == j)
        halfThickness = opticalThickness[i] / 2;
        heatFlux = cellIntensity[i] * (1 - 2 * exponential_integral(halfThickness));

        // Walk away from cell i (along negative coordinate)
        opticalDistance = halfThickness;
        for (auto j = long(i) - 1; j >= 0; --j) {
            heatFlux += cellIntensity[j] * (
                exponential_integral(
                    opticalDistance) -
                    exponential_integral(opticalDistance + opticalThickness[j])
            );
            opticalDistance += opticalThickness[j];
        }
        // Walk away from cell i (along positive coordinate)
        opticalDistance = halfThickness;
        for (auto j = i + 1; j < numCells; ++j) {
            heatFlux += cellIntensity[j] * (
                exponential_integral(
                    opticalDistance) -
                    exponential_integral(opticalDistance + opticalThickness[j])
            );
            opticalDistance += opticalThickness[j];
        }

        heating[i] = heatFlux;
    }

    return heating;
}
