// Approximate element size
lc = 0.0001;
zmax = 0.01;  // Extrusion length

Mesh.RecombinationAlgorithm = 1;
Mesh.RecombineAll = 1;

// Points that contour the copper support
Point(1) = {0.0, 0.0, 0.0, lc};
Point(2) = {0.00485, 0.0, 0.0, lc};
Point(3) = {0.00485, 0.00475, 0.0, lc}; // first point of circle
Point(4) = {0.0096, 0.00475, 0.0, lc}; // center of circle
Point(5) = {0.01435, 0.00475, 0.0, lc}; // second point of circle
Point(6) = {0.01435, 0.0, 0.0, lc};
Point(7) = {0.0192, 0.0, 0.0, lc};
Point(8) = {0.0192, 0.023, 0.0, lc};
Point(9) = {0.0172, 0.023, 0.0, lc}; // this point also contours the stacked HTS tape
Point(10) = {0.002, 0.023, 0.0, lc}; // this point also contours the stacked HTS tape
Point(11) = {0.0, 0.023, 0.0, lc};

// Points that contour the stacked HTS tape
Point(12) = {0.002, 0.011, 0.0, lc};
Point(13) = {0.0172, 0.011, 0.0, lc};

// Center of HTS tape (for air region)
Point(14) = {0.0096, 0.017, 0.0, lc};

// Air region: cylindrical region with radius 45mm centered on HTS tape
Point(15) = {0.0096 + 0.045, 0.017, 0.0, lc};   // right
Point(16) = {0.0096, 0.017 + 0.045, 0.0, lc};   // top
Point(17) = {0.0096 - 0.045, 0.017, 0.0, lc};   // left
Point(18) = {0.0096, 0.017 - 0.045, 0.0, lc};   // bottom

// Lines for copper support contour
Line(1) = {1, 2};      // bottom edge
Line(2) = {2, 3};      // vertical to circle start
Circle(3) = {5, 4, 3}; // circular arc (center at point 4)
Line(4) = {5, 6};      // vertical from circle end
Line(5) = {6, 7};      // bottom right horizontal
Line(6) = {7, 8};      // right vertical
Line(7) = {8, 9};      // top right horizontal
Line(8) = {9, 10};     // top edge (shared with HTS tape)
Line(9) = {10, 11};    // top left horizontal
Line(10) = {11, 1};    // left vertical

// Lines for stacked HTS tape (rectangle inside copper support)
// HTS tape uses points 9, 10 (top), 12, 13 (bottom)
Line(11) = {13, 9};    // right edge of HTS tape
Line(12) = {12, 13};   // bottom edge of HTS tape
Line(13) = {10, 12};   // left edge of HTS tape

// Air region: outer circle (4 quadrants)
Circle(21) = {15, 14, 16};  // quadrant 1: right to top
Circle(22) = {16, 14, 17};  // quadrant 2: top to left
Circle(23) = {17, 14, 18};  // quadrant 3: left to bottom
Circle(24) = {18, 14, 15};  // quadrant 4: bottom to right

// Curve Loops
// Copper support outer boundary (counter-clockwise)
Curve Loop(101) = {1, 2, -3, 4, 5, 6, 7, -11, -12, -13, 9, 10};

// Stacked HTS tape rectangle (counter-clockwise)
Curve Loop(102) = {13, 12, 11, 8};

// Air region: outer circle minus copper support boundary (counter-clockwise)
Curve Loop(103) = {21, 22, 23, 24, -1, -2, 3, -4, -5, -6, -7, -8, -9, -10};

// Plane Surfaces
// Copper support surface with HTS tape as hole
Plane Surface(201) = {101};

// Stacked HTS tape surface
Plane Surface(202) = {102};

// Air region surface
Plane Surface(203) = {103};

// Transfinite meshing for structured grid - 1mm element size
// Line 1: 4.85mm -> 5 elements -> 6 nodes
// Line 2: 4.85mm -> 5 elements -> 6 nodes
// Circle 3: semicircle r=4.75mm, length=14.92mm -> 15 elements -> 16 nodes
// Line 4: 4.85mm -> 5 elements -> 6 nodes
// Line 5: 4.85mm -> 5 elements -> 6 nodes
// Line 6: 23mm -> 23 elements -> 24 nodes
// Line 7: 2mm -> 2 elements -> 3 nodes
// Line 8: 15.2mm -> 15 elements -> 16 nodes
// Line 9: 2mm -> 2 elements -> 3 nodes
// Line 10: 23mm -> 23 elements -> 24 nodes
// Line 11: 12mm -> 12 elements -> 13 nodes
// Line 12: 15.2mm -> 15 elements -> 16 nodes
// Line 13: 12mm -> 12 elements -> 13 nodes
Transfinite Line {1, 2, 4, 5} = 6;
Transfinite Line {6, 10} = 24;
Transfinite Line {7, 9} = 3;
Transfinite Line {8} = 16;
Transfinite Line {11, 13} = 13;
Transfinite Line {12} = 16;
Transfinite Curve {3} = 16;  // copper support circular arc
// Air region circle arcs: quarter circles r=45mm, length=70.69mm -> 71 elements -> 72 nodes
Transfinite Curve {21, 22, 23, 24} = 72;

Transfinite Surface {202}; // too many corners cannot transfinite 201 and 203
Recombine Surface {201, 202, 203};

// Extrusion along z-axis
Extrude {0, 0, zmax} { Surface{201}; Layers{10}; Recombine; }
Extrude {0, 0, zmax} { Surface{202}; Layers{10}; Recombine; }
Extrude {0, 0, zmax} { Surface{203}; Layers{10}; Recombine; }

// Physical groups for volumes
Physical Volume("copper_support") = {1};
Physical Volume("hts_tape") = {2};
Physical Volume("air") = {3};

// Physical groups for surfaces (inlet/outlet)
// Copper support
Physical Surface("copper_support_inlet") = {201};
Physical Surface("copper_support_outlet") = {265}; // extruded top surface

// Stacked HTS tape
Physical Surface("hts_tape_inlet") = {202};
Physical Surface("hts_tape_outlet") = {287};

// Air region
Physical Surface("air_inlet") = {203};
Physical Surface("air_outlet") = {359};
Physical Surface("air_external") = {306,310,314,318};

// Physical groups for side surfaces (extruded sides)
Physical Surface("copper_support_sides") = {220,224,228,232,236,240,244,260,264};
Physical Surface("hts_tape_sides") = {248,252,256,286};