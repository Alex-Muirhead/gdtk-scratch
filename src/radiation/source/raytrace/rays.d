module rays;

import std.format;
import std.math;

import geom.elements.vector3 : Vector3, wedge2D, dot;
import nm.number : number;

interface Ray {
    void walkForward(number length);
    bool intersect(Vector3 vertexOne, Vector3 vertexTwo, out number length);
    Vector3 currentPoint();
}

class CartesianRay : Ray {
public:
    Vector3 terminus;
    Vector3 tangent;
    number rayCoord;

    this(Vector3 tangent, Vector3 terminus) {
        this.tangent = tangent;
        this.terminus = terminus;
        this.rayCoord = 0.0;
    }

    /**
     * Moves the ray forward by the specified length.
     *
     * Params:
     *   length = The distance to move forward (number).
     */
    void walkForward(number length) {
        rayCoord += length;
    }

    Vector3 currentPoint() {
        return terminus + rayCoord * tangent;
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
        Vector3 toVertex = this.currentPoint() - vertexOne;

        number alignment = wedge2D(faceTangent, tangent);
        number s = wedge2D(toVertex, faceTangent) / alignment;
        number t = wedge2D(toVertex, tangent) / alignment;

        debug {
            //         writeln(format("[TRACE] t = %.5f, s = %.5f, alignment = %.2f", t, s, alignment));
        }

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
    number rayCoord;
    number semiMajor; // a
    number semiMinor; // b
    number linearEccentricity; // c

private:
    Vector3 positiveAsymptote;
    Vector3 negativeAsymptote;

public:
    this(number semiMajor, number semiMinor, Vector3 center, number startCoord) {
        this.center = center;
        this.rayCoord = startCoord;

        this.semiMajor = semiMajor;
        this.semiMinor = semiMinor;
        this.linearEccentricity = sqrt(semiMajor ^^ 2 + semiMinor ^^ 2);
        this.positiveAsymptote = Vector3(+semiMinor, semiMajor);
        this.negativeAsymptote = Vector3(-semiMinor, semiMajor);
    }

    /**
     * Constructs a HyperbolicRay instance.
     *
     * Params:
     *   tangent = The tangent vector (Vector3).
     *   terminus = The origin point of the ray (Vector3).
     */
    this(Vector3 tangent, Vector3 terminus) {
        number eccentricity = 1 / sqrt(tangent.y ^^ 2 + tangent.z ^^ 2);
        assert(eccentricity >= 1.0, format("\nEccentricity = %+.3e < 1.0. Not a hyperbola.", eccentricity));

        semiMajor = fabs(terminus.y * tangent.z - terminus.z * tangent.y) * eccentricity; // = a
        semiMinor = tangent.x * semiMajor * eccentricity; // = b
        linearEccentricity = semiMajor * eccentricity; // = c

        rayCoord = (terminus.y * tangent.y + terminus.z * tangent.z) / (tangent.y ^^ 2 + tangent.z ^^ 2);

        Vector3 saddle = terminus - tangent * rayCoord;

        center = Vector3(saddle.x, 0.0); // = f_0

        // Positive and Negative are defined relative to tangent vector.
        // The arc-coordinate is defined in the same frame.
        positiveAsymptote = Vector3(+semiMinor, semiMajor);
        negativeAsymptote = Vector3(-semiMinor, semiMajor);
    }

    override string toString() const {
        return format("HyperbolicRay(a: %.2g, b: %.2g, f0: %s)", semiMajor, semiMinor, center);
    }

    /** 
     * Evaluate the tangent to the hyperbola. 
     *
     * Params:
     *   rayCoord = The rayCoord to evaluate the tangent at (number). This is s
     * Returns: 
     */
    Vector3 localTangent(number rayCoord) {
        number normCoord = rayCoord / linearEccentricity; // This is s'
        return Vector3(semiMinor, semiMajor * normCoord / sqrt(1 + normCoord ^^ 2));
    }

    Vector3 currentPoint() {
        number normCoord = rayCoord / linearEccentricity; // s'
        return center + Vector3(semiMinor * normCoord, semiMajor * sqrt(1 + normCoord ^^ 2));
    }

    /**
     * Moves the ray forward by the specified arc length.
     *
     * Params:
     *   arcLength = The distance to move forward (number).
     */
    void walkForward(number arcLength) {
        rayCoord += arcLength;
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
    bool intersect(Vector3 vertexOne, Vector3 vertexTwo, out number length) {
        // Calculate the face tangent vector
        Vector3 faceTangent = vertexTwo - vertexOne;
        // Calculate the vector from the center to the first vertex
        Vector3 toVertex = vertexOne - center;

        number positiveBeta = wedge2D(toVertex, faceTangent) / wedge2D(positiveAsymptote, faceTangent);
        number negativeBeta = wedge2D(toVertex, faceTangent) / wedge2D(negativeAsymptote, faceTangent);

        number invAreaRatio = 1 / (positiveBeta * negativeBeta);
        // Needs an area ratio greater than 1 to have any intersection
        if (invAreaRatio > 1) {
            return false;
        }
        // Potentially, if areaRatio >> 1, we can drop back to linear?

        number bias = sqrt(1 - invAreaRatio);

        // Try both solutions using static foreach
        // NOTE: `static foreach` has the limitation that you can't use `continue` etc
        // OH WAIT! We only use the `sign` in the two exponential solutions, but we can
        // Throw away one of them entirely
        // ...but there might be two solutions still if sqrt(determinant) < findAName
        foreach (sign; [+1, -1]) {
            {
                // Calculate the exponential, arc, Cartesian, and linear coordinates at the intercept
                number expHypAngle = positiveBeta * (1 + sign * bias);
                if (expHypAngle < 0) {
                    continue;
                }

                number interceptCoord = (expHypAngle ^^ 2 - 1) / (2 * expHypAngle) * linearEccentricity;
                number crossingDirection = wedge2D(localTangent(interceptCoord), faceTangent);
                if (crossingDirection < 0) {
                    continue;
                }

                Vector3 cartesian = center + positiveAsymptote * (expHypAngle / 2) + negativeAsymptote / (
                    expHypAngle * 2);
                Vector3 translation = cartesian - vertexOne;
                number linear = dot(translation, faceTangent) / dot(faceTangent, faceTangent);

                // Check if the intersection point is within the line segment
                // and ahead of the current ray position
                if (linear >= 0 && linear < 1 && interceptCoord > rayCoord) {
                    length = interceptCoord - rayCoord;
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

@("Hyperbolic Intercept - Near-Horizontal Edge")
unittest {
    import fluent.asserts;

    number semiMajor = 0.5;
    number semiMinor = 0.5;
    Vector3 center = Vector3(0.0, 0.0);

    Ray ray = new HyperbolicRay(semiMajor, semiMinor, center, 0.0);

    // Define edge vertices
    Vector3 vertex0 = Vector3(1.2, 0.8);
    Vector3 vertex1 = Vector3(0.5, 0.9);

    number arcLength;

    // The crossing should be "positive" (aka exiting the cell)
    bool success = ray.intersect(vertex0, vertex1, arcLength);
    Assert.equal(success, true);

    // Flip the edge order, and the cross should prevent a solution
    success = ray.intersect(vertex1, vertex0, arcLength);
    Assert.equal(success, false);
}
