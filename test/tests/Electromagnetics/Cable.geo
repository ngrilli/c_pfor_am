// Approximate element size
lc = 0.0001;
zmax = 0.01;  // Extrusion length

Mesh.RecombinationAlgorithm = 1;
Mesh.RecombineAll = 1;

// Points that contour the copper support
Point(1) = {0.0, 0.0, 0.0, lc};
Point(2) = {0.00485, 0.0, 0.0, lc};
Point(3) = {0.00485, 0.00485, 0.0, lc}; // first point of circle
Point(4) = {0.0096, 0.00485, 0.0, lc}; // center of circle
Point(5) = {0.01435, 0.00485, 0.0, lc}; // second point of circle
Point(6) = {0.01435, 0.0, 0.0, lc};
Point(7) = {0.0192, 0.0, 0.0, lc};
Point(8) = {0.0192, 0.023, 0.0, lc};
Point(9) = {0.0172, 0.023, 0.0, lc}; // this point also contours the stacked HTS tape
Point(10) = {0.002, 0.023, 0.0, lc}; // this point also contours the stacked HTS tape
Point(11) = {0.0, 0.023, 0.0, lc};

// Points that contour the stacked HTS tape
Point(12) = {0.002, 0.011, 0.0, lc};
Point(13) = {0.0172, 0.011, 0.0, lc};

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

// Curve Loops
// Copper support outer boundary (counter-clockwise)
Curve Loop(101) = {1, 2, -3, 4, 5, 6, 7, -11, -12, -13, 9, 10};

// Stacked HTS tape rectangle (counter-clockwise)
Curve Loop(102) = {13, 12, 11, 8};

// Plane Surfaces
// Copper support surface with HTS tape as hole
Plane Surface(201) = {101};

// Stacked HTS tape surface
Plane Surface(202) = {102};

// Transfinite meshing for structured grid
Transfinite Line {1, 4, 5, 7, 9, 10, 11, 13} = 50;
Transfinite Line {2, 6} = 20;
Transfinite Curve {3} = 30;  // circular arc
Transfinite Line {8, 12, 14} = 50;  // top/bottom edges shared

Transfinite Surface {201, 202};
Recombine Surface {201, 202};

// Extrusion along z-axis
Extrude {0, 0, zmax} { Surface{201}; Layers{10}; Recombine; }
Extrude {0, 0, zmax} { Surface{202}; Layers{10}; Recombine; }

// Physical groups for volumes
Physical Volume("copper_support") = {1};
Physical Volume("hts_tape") = {2};

// Physical groups for surfaces (inlet/outlet)
// Copper support
Physical Surface("copper_support_inlet") = {201};
Physical Surface("copper_support_outlet") = {201 + 2000}; // extruded top surface

// Stacked HTS tape
Physical Surface("hts_tape_inlet") = {202};
Physical Surface("hts_tape_outlet") = {202 + 2000};

// Physical groups for side surfaces (extruded sides)
Physical Surface("copper_support_sides") = {2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008, 2009, 2010};
Physical Surface("hts_tape_sides") = {2011, 2012, 2013, 2014};