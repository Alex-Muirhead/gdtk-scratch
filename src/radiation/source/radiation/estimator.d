module radiation.estimator;

import std.math;

import nm.number;
import geom.elements.vector3 : Vector3, wedge2D, dot;

import lmr.fvcell : FVCell;


// NOTE: We use a generic `T` here.
// What I really need is for `T` to be a "Field Element" (in the mathematical sense),
// such that element addition and scalar multiplication (and hence integer exponentiation)
// are correctly defined. Vectors should be applicable within this.
// It also follows there should be an identity element under `T`

struct SimpleEstimator(T) {
    size_t zeroth; // zeroth-moment, or the count
    T first; // first-moment, or the sum
    T second; // second-moment, or the sum of squares

    SimpleEstimator opBinary(string op : "+")(T value) {
        return SimpleEstimator(zeroth + 1, first + value, second + value ^^ 2);
    }

    void opOpAssign(string op : "+")(T value) {
        zeroth += 1;
        first += value;
        second += value^^2;
    }

    @property
    size_t sample_count() {
        return zeroth;
    }

    @property
    T sample_mean() {
        return first / zeroth;
    }

    @property
    T sample_var() {
        return second / zeroth - (first / zeroth) ^^ 2;
    }

    @property
    T value() {
        return sample_mean();
    }

    @property
    T stddev() {
        return sqrt(sample_var() / sample_count());
    }
}

@("SimpleEstimator -- Basic functionality")
unittest {
    import fluent.asserts;
    import std.stdio;

    auto est = SimpleEstimator!number(0, 0.0, 0.0);
    auto updated = est + 2.0;
    updated += 4.0;

    Assert.equal(updated.first, 6.0);

    Assert.equal(updated.sample_count, 2);
    Assert.equal(updated.sample_mean, 3.0);
    Assert.equal(updated.sample_var, 1.0);
}

struct Location {
    size_t cellID;
    Vector3 absolute;
    // Barycentric values?
}

struct Crossing {
    FVCell cell;  // or maybe cellID?
    number length;
}

class RayPath {

}

// void tracing() {
//     // Generate this from cell volume and temperature
//     Location origin = ...;
//     auto packetEnergy = ...;
//     auto ray = ...;
//     auto rayPath = ;
//     foreach (crossing; rayPath) {
//         // We can mut rayPath here too (oops bad ownership?)
//         // Do our Monte-Carlo integration
//         crossing.cell.fs.Qrad += ... ;

//         // Handle boundaries here
//         if (boundary.type == "GhostType") {
//             rayPath.teleport(boundary_stuff);
//         }
//     }
// }
