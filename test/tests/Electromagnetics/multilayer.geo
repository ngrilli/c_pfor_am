lc = 0.2;

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
Point(12) = {-0.003, 0.000047, 0.0, lc};
Point(13) = {0.003, 0.000047, 0.0, lc};

// copper bottom
Line(101) = {2, 3};
Line(102) = {3, 4};
Line(103) = {4, 5};
Line(104) = {5, 2};

// substrate
Line(105) = {4, 6};
Line(106) = {6, 7};
Line(107) = {7, 5};

// inner circle
//Point(2) = {1.0, 0.0, 0.0, lc};
//Point(3) = {0.0, 1.0, 0.0, lc};
//Point(4) = {-1.0, 0.0, 0.0, lc};
//Point(5) = {0.0, -1.0, 0.0, lc};

// outer circle
//Point(6) = {2.0, 0.0, 0.0, lc};
//Point(7) = {0.0, 2.0, 0.0, lc};
//Point(8) = {-2.0, 0.0, 0.0, lc};
//Point(9) = {0.0, -2.0, 0.0, lc};

// inner circle
//Circle(11) = {2, 1, 5};
//Circle(12) = {5, 1, 4};
//Circle(13) = {4, 1, 3};
//Circle(14) = {3, 1, 2};

// outer circle
//Circle(15) = {6, 1, 7};
//Circle(16) = {7, 1, 8};
//Circle(17) = {8, 1, 9};
//Circle(18) = {9, 1, 6};

Curve Loop(201) = {101, 102, 103, 104};
Curve Loop(202) = {-103, 105, 106, 107};

//Curve Loop(21) = {11, 12, 13, 14};
//Curve Loop(22) = {15, 16, 17, 18};
//Curve Loop(23) = {-11, -12, -13, -14};

Plane Surface(201) = {201};
Plane Surface(202) = {202};

//Plane Surface(31) = {21, 22};
//Plane Surface(32) = {23};

Recombine Surface{201};
Recombine Surface{202};

//Recombine Surface{32};

Extrude {0,0,0.010} { Surface{201}; Layers{20}; Recombine; }
Extrude {0,0,0.010} { Surface{202}; Layers{20}; Recombine; }

//Extrude {0,0,5} { Surface{32}; Layers{25}; Recombine; }

//Physical Surface("outer_back") = {31};
//Physical Surface("inner_back") = {32};
//Physical Surface("outer_front") = {74};
//Physical Surface("inner_front") = {96};
//Physical Surface("outer_circle") = {61, 65, 69, 73};

//Physical Volume("inner_block") = {2};
//Physical Volume("outer_block") = {1};