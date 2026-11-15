SetFactory("OpenCASCADE");

lc = 0.000001;

Mesh.Algorithm = 6;
Mesh.RecombinationAlgorithm = 1;
Mesh.RecombineAll = 1;
Mesh.CharacteristicLengthMin = lc;
Mesh.CharacteristicLengthMax = 100*lc;

// Geometry dimensions
// Inner rectangular core: x in [-0.003, 0.003], y in [-0.000046, 0.000047]
// Outer cylinder: radius 0.009
// Extrusion in z: [0, lc]

// Layer boundaries in y (from bottom to top)
y0 = -0.000046;  // copper_bottom bottom
y1 = -0.000026;  // copper_bottom top / substrate bottom
y2 =  0.000024;  // substrate top / tape bottom
y3 =  0.000025;  // tape top / silver bottom
y4 =  0.000027;  // silver top / copper_top bottom
y5 =  0.000047;  // copper_top top

x_min = -0.003;
x_max =  0.003;
z_min =  0.0;
z_max =  5*lc;

R_outer = 0.009;

// Outer cylinder (full height)
Cylinder(1) = {0.0, 0.0, 0.0, 0.0, 0.0, z_max, R_outer, 2*Pi};

// Inner rectangular core boxes for each layer
// copper_bottom
Box(2) = {x_min, y0, z_min, x_max - x_min, y1 - y0, z_max - z_min};
// substrate
Box(3) = {x_min, y1, z_min, x_max - x_min, y2 - y1, z_max - z_min};
// tape
Box(4) = {x_min, y2, z_min, x_max - x_min, y3 - y2, z_max - z_min};
// silver
Box(5) = {x_min, y3, z_min, x_max - x_min, y4 - y3, z_max - z_min};
// copper_top
Box(6) = {x_min, y4, z_min, x_max - x_min, y5 - y4, z_max - z_min};

// BooleanFragments to separate all volumes
// First fragment the outer cylinder with all inner boxes
v() = BooleanFragments{ Volume{1}; Delete; }{ Volume{2, 3, 4, 5, 6}; Delete; };

// After BooleanFragments, volumes are re-indexed. We need to identify them.

// Physical Volumes
Physical Volume("copper_bottom", 1) = {2};
Physical Volume("substrate", 2) = {3};
Physical Volume("tape", 3) = {4};
Physical Volume("silver", 4) = {5};
Physical Volume("copper_top", 5) = {6};
Physical Volume("outer_circle", 6) = {7};

// Physical Surfaces - Inlet (z=0) and Outlet (z=5*lc) faces

// Copper bottom (volume 2): inlet=surface of vol 2 at z=0, outlet=surface of vol 2 at z=lc
Physical Surface("copper_bottom_inlet", 1) = {17};
Physical Surface("copper_bottom_outlet", 2) = {18};

// Substrate (volume 3): inlet/outlet
Physical Surface("substrate_inlet", 3) = {20};
Physical Surface("substrate_outlet", 4) = {21};

// Tape (volume 4): inlet/outlet
Physical Surface("tape_inlet", 5) = {23};
Physical Surface("tape_outlet", 6) = {24};

// Silver (volume 5): inlet/outlet
Physical Surface("silver_inlet", 7) = {26};
Physical Surface("silver_outlet", 8) = {27};

// Copper top (volume 6): inlet/outlet
Physical Surface("copper_top_inlet", 9) = {28};
Physical Surface("copper_top_outlet", 10) = {29};

// Outer circle (volume 1): inlet/outlet (annular faces at z=0 and z=lc)
Physical Surface("outer_circle_inlet", 11) = {3};
Physical Surface("outer_circle_outlet", 12) = {2};

// External surface of outer circle (cylindrical surface)
Physical Surface("outer_circle_external", 13) = {1};

// right side internal surfaces
Physical Surface("copper_bottom_right", 14) = {5};
Physical Surface("substrate_right", 15) = {6};
Physical Surface("tape_right", 16) = {7};
Physical Surface("silver_right", 17) = {8};
Physical Surface("copper_top_right", 18) = {9};

// left side internal surfaces
Physical Surface("copper_bottom_left", 19) = {15};
Physical Surface("substrate_left", 20) = {14};
Physical Surface("tape_left", 21) = {13};
Physical Surface("silver_left", 22) = {12};
Physical Surface("copper_top_left", 23) = {11};

// bottom/top side internal surfaces
Physical Surface("copper_bottom_bottom", 24) = {4};
Physical Surface("substrate_bottom", 25) = {16};
Physical Surface("tape_bottom", 26) = {19};
Physical Surface("silver_bottom", 27) = {22};
Physical Surface("copper_top_bottom", 28) = {25};
Physical Surface("copper_top_top", 29) = {10};

// Mesh settings
Mesh.MshFileVersion = 2.2;
Mesh.Binary = 0;
Mesh.RecombineAll = 1;

// Transfinite settings for structured mesh (similar to original)
// Note: Transfinite doesn't work directly with OpenCASCADE BooleanFragments
// The mesh will be unstructured but we can set characteristic length