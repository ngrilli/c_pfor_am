lc = 0.2;

Point(1) = {0.0, 0.0, 0.0, lc}; // center

// inner circle (radius 1.0)
Point(2) = {1.0, 0.0, 0.0, lc};
Point(3) = {0.0, 1.0, 0.0, lc};
Point(4) = {-1.0, 0.0, 0.0, lc};
Point(5) = {0.0, -1.0, 0.0, lc};

// outer circle (radius 2.0)
Point(6) = {2.0, 0.0, 0.0, lc};
Point(7) = {0.0, 2.0, 0.0, lc};
Point(8) = {-2.0, 0.0, 0.0, lc};
Point(9) = {0.0, -2.0, 0.0, lc};

// inner circle
Circle(11) = {2, 1, 5};
Circle(12) = {5, 1, 4};
Circle(13) = {4, 1, 3};
Circle(14) = {3, 1, 2};

// outer circle
Circle(15) = {6, 1, 7};
Circle(16) = {7, 1, 8};
Circle(17) = {8, 1, 9};
Circle(18) = {9, 1, 6};

Curve Loop(21) = {11, 12, 13, 14};
Curve Loop(22) = {15, 16, 17, 18};

Plane Surface(31) = {21, 22};

Recombine Surface{31};

// Physical boundaries (2D: Physical Curve for boundaries, Physical Surface for domain)
Physical Curve("inner_circle") = {11, 12, 13, 14};
Physical Curve("outer_circle") = {15, 16, 17, 18};
Physical Surface("domain") = {31};