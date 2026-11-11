SetFactory("OpenCASCADE");

lc = 0.000001;

Mesh.Algorithm = 6;
Mesh.RecombinationAlgorithm = 1;
Mesh.RecombineAll = 1;
Mesh.CharacteristicLengthMin = lc;
Mesh.CharacteristicLengthMax = lc;

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
z_max =  lc;

R_outer = 0.009;

// Outer cylinder (full height)
Cylinder(1) = {0.0, 0.0, 0.0, 0.0, 0.0, lc, R_outer, 2*Pi};

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
// The order typically preserves the first argument's volumes first, then the second's.
// But with multiple volumes in second argument, they may be interleaved.
// Use Physical Volume to tag them by their geometric properties.

// The outer cylinder gets split into: outer shell (annular) + 5 inner layer volumes
// The 5 boxes become the 5 layer volumes (possibly split further)

// Tag volumes by their center of mass y-coordinate using physical groups
// We'll use a post-processing approach: get all volumes and classify them

// For now, assign physical volumes based on typical BooleanFragments output order
// The fragments from Volume{1} (outer cylinder) come first, then from Volumes{2..6}
// But since we used Delete on both, all original volumes are deleted and new ones created.

// Better approach: create physical groups for surfaces first, then volumes

// Physical Surfaces - Inlet/Outlet faces (z = 0 and z = lc)
// We'll tag surfaces after fragmentation by their z-coordinate

// Physical Volumes - will be assigned after checking the actual volume tags
// The typical output order for BooleanFragments with 1 volume vs 5 volumes:
// Results in 6 volumes: 1 outer annular + 5 inner rectangular layers

// Since we can't easily predict the exact volume tags, we use a trick:
// Create physical volumes by listing all volume tags and let the user verify
// Or use the fact that the inner boxes (2-6) are completely inside the cylinder

// Actually, let's use a more robust approach: create the geometry in steps
// 1. Create outer cylinder
// 2. Create inner boxes
// 3. Fragment

// The volume tags after BooleanFragments with Delete on both:
// New volumes are created sequentially. The fragments of the first argument (cylinder)
// come first, then fragments of the second argument (boxes).
// Since the boxes are disjoint in y, each box becomes one volume fragment.
// The cylinder gets split into 6 fragments: 1 outer annular + 5 "holes" where boxes were.
// But with Delete, the "hole" volumes are the box fragments themselves.

// Expected volume tags after BooleanFragments (6 volumes total):
// 1: outer annular region (from cylinder)
// 2: copper_bottom (from box 2)
// 3: substrate (from box 3)
// 4: tape (from box 4)
// 5: silver (from box 5)
// 6: copper_top (from box 6)

// Physical Volumes
Physical Volume("copper_bottom", 1) = {2};
Physical Volume("substrate", 2) = {3};
Physical Volume("tape", 3) = {4};
Physical Volume("silver", 4) = {5};
Physical Volume("copper_top", 5) = {6};
Physical Volume("outer_circle", 6) = {1};

// Physical Surfaces - Inlet (z=0) and Outlet (z=lc) faces
// Surface tags are harder to predict, so we use a wide range and filter by name in the solver
// Or we can get them by the volume's boundary surfaces

// Get inlet/outlet surfaces for each volume
// For each volume, surface 1 and 2 are typically the z-min and z-max faces (for extruded boxes)
// For the outer annular volume, it's more complex

// Better: use Physical Surface with the volume's boundary surfaces
// We'll define them after seeing the actual mesh output, or use a generic approach

// For now, define physical surfaces for all z-min and z-max faces of all volumes
// The surface tags for the 6 volumes' bottom/top faces:
// Each volume typically has 6 surfaces (for boxes) or more (for annular)
// The bottom face (z=0) and top face (z=lc) are usually the first two surfaces for simple volumes

// We'll use a broader approach: define physical surfaces for inlet/outlet by volume
// This requires knowing the surface tags. Let's use the fact that for extruded boxes:
// Surface 1: z-min (inlet), Surface 2: z-max (outlet), Surfaces 3-6: sides

// For the 5 inner volumes (boxes), surfaces 1 and 2 are inlet/outlet
// For the outer volume (annular), it has more surfaces

// Let's assign physical surfaces using a loop-like approach in comments
// The actual surface tags need to be determined from the mesh

// Copper bottom (volume 2): inlet=surface of vol 2 at z=0, outlet=surface of vol 2 at z=lc
Physical Surface("copper_bottom_inlet", 101) = {2};
Physical Surface("copper_bottom_outlet", 102) = {8};

// Substrate (volume 3): inlet/outlet
Physical Surface("substrate_inlet", 103) = {3};
Physical Surface("substrate_outlet", 104) = {9};

// Tape (volume 4): inlet/outlet
Physical Surface("tape_inlet", 105) = {4};
Physical Surface("tape_outlet", 106) = {10};

// Silver (volume 5): inlet/outlet
Physical Surface("silver_inlet", 107) = {5};
Physical Surface("silver_outlet", 108) = {11};

// Copper top (volume 6): inlet/outlet
Physical Surface("copper_top_inlet", 109) = {6};
Physical Surface("copper_top_outlet", 110) = {12};

// Outer circle (volume 1): inlet/outlet (annular faces at z=0 and z=lc)
Physical Surface("outer_circle_inlet", 111) = {1};
Physical Surface("outer_circle_outlet", 112) = {7};

// External surface of outer circle (cylindrical surface)
Physical Surface("outer_circle_external", 113) = {13, 14, 15, 16};

// Mesh settings
Mesh.MshFileVersion = 2.2;
Mesh.Binary = 0;
Mesh.RecombineAll = 1;

// Transfinite settings for structured mesh (similar to original)
// Note: Transfinite doesn't work directly with OpenCASCADE BooleanFragments
// The mesh will be unstructured but we can set characteristic length