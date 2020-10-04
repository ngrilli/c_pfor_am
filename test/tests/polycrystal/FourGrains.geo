// Nicolo Grilli
// National University of Singapore
// 2 October 2020

// Square geometry, 4 grains

Point(1) = {0, 0, 0, 0.25};
Point(2) = {0.5, 0, 0, 0.25};
Point(3) = {1, 0, 0, 0.25};
Point(4) = {0, 0.5, 0, 0.25};
Point(5) = {0.5, 0.5, 0, 0.25};
Point(6) = {1, 0.5, 0, 0.25};
Point(7) = {0, 1, 0, 0.25};
Point(8) = {0.5, 1, 0, 0.25};
Point(9) = {1, 1, 0, 0.25};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 6};
Line(4) = {6, 9};
Line(5) = {9, 8};
Line(6) = {8, 7};
Line(7) = {7, 4};
Line(8) = {4, 1};
Line(9) = {5, 2};
Line(10) = {6, 5};
Line(11) = {8, 5};
Line(12) = {5, 4};

// Line loops

Line Loop(1) = {1, -9, 12, 8};
Line Loop(2) = {2, 3, 10, 9};
Line Loop(3) = {4, 5, 11, -10};
Line Loop(4) = {6, 7, -12, -11};

Plane Surface(11) = {1};
Plane Surface(12) = {2};
Plane Surface(13) = {3};
Plane Surface(14) = {4};

// Transfinite Line

Transfinite Line {1,2,3,4,5,6,7,8,9,10,11,12} = 3 Using Progression 1;

Transfinite Surface {11};
Transfinite Surface {12};
Transfinite Surface {13};
Transfinite Surface {14};

// Recombine Surface

Recombine Surface {11};
Recombine Surface {12};
Recombine Surface {13};
Recombine Surface {14};

Extrude {0, 0, 0.25} {
  Surface{11,12,13,14};
  Layers{1};
  Recombine;
}

Physical Volume("grain1") = {1};
Physical Volume("grain2") = {2};
Physical Volume("grain3") = {3};
Physical Volume("grain4") = {4};

Physical Surface("Back") = {13, 14, 11, 12};
Physical Surface("Bottom") = {23, 45};
Physical Surface("Left") = {35, 93};






