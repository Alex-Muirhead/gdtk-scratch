module radiation.ray.directions;

import std.format;
import std.math;

import geom.elements.vector3 : Vector3, wedge2D, dot;
import nm.number : number;

// External packages

import mir.random.engine : Random;
import mirRandom = mir.random.engine;
import mir.random.ndvariable : sphereVar;
import mir.random.variable : uniformVar;

class Stochastic {
public:
    int dimensions;

    this(int dimensions) {
        this.dimensions = dimensions;
    }
}

class Deterministic {
private:
    static PHI = (1.0 + sqrt(5.0)) / 2.0;
    static OFFSET = PI / 4;
    int i;
    int N;
public:
    uint dimensions;
    uint samples;

    this(uint dimensions, uint samples) {
        assert(samples % 2 == 1, "Number of samples must be odd.");
        this.dimensions = dimensions;
        this.samples = samples;
        this.N = (samples - 1) / 2;
        this.i = -N;
    }

    // foreach impl block
    bool empty() const {
        return i == N+1;
    }

    void popFront() {
        ++i;
    }

    Vector3 front() const {
        number lon = 2*PI*i / PHI + OFFSET; // Add offset to avoid collision with YZ-plane
        number z = 2.0*i / samples; // lat = arcsin(z);
        number comp = sqrt(1 - z*z);
        return Vector3(x: comp*sin(lon), y: comp*cos(lon), z: z);
    }
}

@("Fibonacci Sphere - PHI Value")
unittest {
    import fluent.asserts;

    Assert.approximately(Deterministic.PHI, 1.618, 1e-3);
}

@("Fibonacci Sphere - Single Point")
unittest {
    import fluent.asserts;
    import std.stdio;

    auto samples = new Deterministic(3, 1);
    Vector3 direction = samples.front();

    Assert.equal(direction.z, 0.0);
    Assert.equal(direction.x, direction.y);

    samples.popFront();
    Assert.equal(samples.empty(), true);
}

@("Fibonacci Sphere - Collect into array")
unittest {
    import fluent.asserts;
    import std.array;
    import std.stdio;

    auto samples = new Deterministic(3, 5);

    Vector3[] directions = array(samples);
    Assert.equal(directions.length, 5);
}
