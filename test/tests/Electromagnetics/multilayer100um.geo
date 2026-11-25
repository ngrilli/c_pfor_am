// 100 micron wide and long elements in the tape

lc = 0.000100;
zmax = 0.000500;

Mesh.RecombinationAlgorithm = 1;
Mesh.RecombineAll = 1;

// center
Point(1) = {0.0, 0.0, 0.0, lc}; 

// copper bottom
Point(2) = {-0.003, -0.000046, 0.0, lc};
Point(3) = {0.003, -0.000046, 0.0, lc};
Point(4) = {0.003, -0.000026, 0.0, lc};
Point(5) = {-0.003, -0.000026, 0.0, lc};

// substrate
Point(6) = {0.003, 0.000024, 0.0, lc};
Point(7) = {-0.003, 0.000024, 0.0, lc};

// tape
Point(8) = {0.003, 0.000025, 0.0, lc};
Point(9) = {-0.003, 0.000025, 0.0, lc};

// silver
Point(10) = {0.003, 0.000027, 0.0, lc};
Point(11) = {-0.003, 0.000027, 0.0, lc};

// copper top
Point(12) = {0.003, 0.000047, 0.0, lc};
Point(13) = {-0.003, 0.000047, 0.0, lc};

// outer circle
Point(14) = {0.009, 0.0, 0.0, lc};
Point(15) = {0.0, 0.009, 0.0, lc};
Point(16) = {-0.009, 0.0, 0.0, lc};
Point(17) = {0.0, -0.009, 0.0, lc};

// copper bottom
Line(101) = {2, 3};
Line(102) = {3, 4};
Line(103) = {4, 5};
Line(104) = {5, 2};

// substrate
Line(105) = {4, 6};
Line(106) = {6, 7};
Line(107) = {7, 5};

// tape
Line(108) = {6, 8};
Line(109) = {8, 9};
Line(110) = {9, 7};

// silver
Line(111) = {8, 10};
Line(112) = {10, 11};
Line(113) = {11, 9};

// copper top
Line(114) = {10, 12};
Line(115) = {12, 13};
Line(116) = {13, 11};

// outer circle
Circle(151) = {14, 1, 15};
Circle(152) = {15, 1, 16};
Circle(153) = {16, 1, 17};
Circle(154) = {17, 1, 14};

Curve Loop(201) = {101, 102, 103, 104};
Curve Loop(202) = {-103, 105, 106, 107};
Curve Loop(203) = {-106, 108, 109, 110};
Curve Loop(204) = {-109, 111, 112, 113};
Curve Loop(205) = {-112, 114, 115, 116};

// outer circle
Curve Loop(206) = {151, 152, 153, 154, -101, -102, -105, -108, -111, -114, -115, -116, -113, -110, -107, -104};

Transfinite Line {101, 103, 106, 109, 112, 115} = 61;
Transfinite Line {102, 104, 114, 116} = 4;
Transfinite Line {105, 107} = 7;
Transfinite Line {108, 110, 111, 113} = 2;

// Number of nodes on each outer circular arc
Transfinite Curve {152, 153} = 21;
Transfinite Curve {151, 154} = 22;

Plane Surface(201) = {201};
Plane Surface(202) = {202};
Plane Surface(203) = {203};
Plane Surface(204) = {204};
Plane Surface(205) = {205};
Plane Surface(206) = {206};

Transfinite Surface {201, 202, 203, 204, 205};

Recombine Surface{201};
Recombine Surface{202};
Recombine Surface{203};
Recombine Surface{204};
Recombine Surface{205};
Recombine Surface{206};

Extrude {0,0,zmax} { Surface{201}; Layers{5}; Recombine; }
Extrude {0,0,zmax} { Surface{202}; Layers{5}; Recombine; }
Extrude {0,0,zmax} { Surface{203}; Layers{5}; Recombine; }
Extrude {0,0,zmax} { Surface{204}; Layers{5}; Recombine; }
Extrude {0,0,zmax} { Surface{205}; Layers{5}; Recombine; }
Extrude {0,0,zmax} { Surface{206}; Layers{5}; Recombine; }

Physical Surface("copper_bottom_inlet") = {201};
Physical Surface("copper_bottom_outlet") = {228};
Physical Surface("substrate_inlet") = {202};
Physical Surface("substrate_outlet") = {250};
Physical Surface("tape_inlet") = {203};
Physical Surface("tape_outlet") = {272};
Physical Surface("silver_inlet") = {204};
Physical Surface("silver_outlet") = {294};
Physical Surface("copper_top_inlet") = {205};
Physical Surface("copper_top_outlet") = {316};
Physical Surface("outer_circle_inlet") = {206};
Physical Surface("outer_circle_outlet") = {398};
Physical Surface("outer_circle_external") = {337, 341, 345, 349};

Physical Surface("copper_bottom_right") = {219};
Physical Surface("substrate_right") = {241};
Physical Surface("tape_right") = {263};
Physical Surface("silver_right") = {285};
Physical Surface("copper_top_right") = {307};
Physical Surface("copper_bottom_left") = {227};
Physical Surface("substrate_left") = {249};
Physical Surface("tape_left") = {271};
Physical Surface("silver_left") = {293};
Physical Surface("copper_top_left") = {315};

Physical Surface("copper_bottom_bottom") = {215};
Physical Surface("substrate_bottom") = {223};
Physical Surface("tape_bottom") = {245};
Physical Surface("silver_bottom") = {267};
Physical Surface("copper_top_bottom") = {289};
Physical Surface("copper_top_top") = {311};

Physical Volume("copper_bottom") = {1};
Physical Volume("substrate") = {2};
Physical Volume("tape") = {3};
Physical Volume("silver") = {4};
Physical Volume("copper_top") = {5};
Physical Volume("outer_circle") = {6};