SetFactory("OpenCASCADE");

lc = 0.1;

Mesh.Algorithm = 6;
Mesh.CharacteristicLengthMin = lc;
Mesh.CharacteristicLengthMax = lc;

Rinner = 1.0;
Router = 2.0;

Cylinder(1) = {0.0, 0.0, 0.0, 0.0, 0.0, 1.0, Rinner, 2*Pi};
Cylinder(2) = {0.0, 0.0, 0.0, 0.0, 0.0, 1.0, Router, 2*Pi};

v() = BooleanFragments{ Volume{1}; Delete; }{ Volume{2}; Delete; };

Physical Surface("inner_circle", 1) = {1};
Physical Surface("outer_circle", 4) = {4};
Physical Surface("inner_back", 3) = {3};
Physical Surface("inner_front", 2) = {2};
Physical Surface("outer_back", 6) = {6};
Physical Surface("outer_front", 5) = {5};

Physical Volume("inner_block", 1) = {1};
Physical Volume("outer_block", 2) = {2};

Mesh.MshFileVersion = 2.2;
Mesh.Binary = 0;