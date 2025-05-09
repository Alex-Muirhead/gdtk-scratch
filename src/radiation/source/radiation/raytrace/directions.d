module directions;

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
public:
    int dimensions;

    this(int dimensions) {
        this.dimensions = dimensions;
    }
}
