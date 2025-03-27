module simple_grid;

import std.getopt;
import std.stdio;

import geom;

void main(string[] args) {

    string workingDir = ".";
    getopt(args, std.getopt.config.stopOnFirstNonOption, "d|dir", &workingDir);

    auto width = 1.0;
    auto height = 10.0;

    // EXAMPLE: Constructing grid in D
    // d -- c
    // |    |
    // a -- b
    auto a = Vector3(0, 0, 0);
    auto b = Vector3(width, 0, 0);
    auto c = Vector3(width, height, 0);
    auto d = Vector3(0, height, 0);
    auto patch = new CoonsPatch(a, b, c, d);

    auto sgrid = new StructuredGrid(patch, 4 + 1, 4 + 1, []);
    writeln("Structured Grid is ", sgrid);

    auto usgrid = new UnstructuredGrid(sgrid, new_label:
        "simple");
    writeln("Unstructured Grid is ", usgrid);
}
