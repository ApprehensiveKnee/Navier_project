// Gmsh project created on Fri Aug  4 23:35:20 2023

SetFactory("OpenCASCADE");

lcar1 = 0.25;
lcar2 = .05;
lcar3 = .025;

// Genarte a box with dimensions
// x = 3.5
// y = 2
// z = 1.5
// left corner at (0,0,0)
Point(1) = {0,0,0, lcar1};
Point(2) = {3.5,0,0, lcar1};
Point(3) = {3.5,2,0, lcar1};
Point(4) = {0,2,0, lcar1};
Point(5) = {0,0,1, lcar1};
Point(6) = {3.5,0,1, lcar1};
Point(7) = {3.5,2,1, lcar1};
Point(8) = {0,2,1, lcar1};

// Now generate the lines
Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};
Line(5) = {5,6};
Line(6) = {6,7};
Line(7) = {7,8};
Line(8) = {8,5};
Line(9) = {1,5};
Line(10) = {2,6};
Line(11) = {3,7};
Line(12) = {4,8};

// Now generate the surfaces
Line Loop(1) = {1,2,3,4};
Plane Surface(1) = {1};
Line Loop(2) = {5,6,7,8};
Plane Surface(2) = {2};
Line Loop(3) = {1,10,-5,-9};
Plane Surface(3) = {3};
Line Loop(4) = {2,11,-6,-10};
Plane Surface(4) = {4};
Line Loop(5) = {3,12,-7,-11};
Plane Surface(5) = {5};
Line Loop(6) = {4,9,-8,-12};
Plane Surface(6) = {6};

// Now generate the surface loop
Surface Loop(1) = {1,2,3,4,5,6};
Volume(1) = {1};


// Generate the cylindrical hole:
// radius 0.2
// height 1
// axis along y
// bottom center at (0.75,0.5,0.5)
Point(9) = {0.75,0.5,0.5,lcar3};
r  = 0.2;

Point(10) = {0.75 + r,0.5,0.5,lcar3};
Point(11) = {0.75,0.5,0.5 + r,lcar3};
Point(12) = {0.75 - r,0.5,0.5,lcar3};
Point(13) = {0.75,0.5,0.5 - r,lcar3};

Circle(13) = {10,9,11};
Circle(14) = {11,9,12};
Circle(15) = {12,9,13};
Circle(16) = {13,9,10};

Line Loop(7) = {13,14,15,16};
Plane Surface(7) = {7};

Extrude {0,1,0} {
  Surface{7};
}
// Assign a mesh size to all the points of all the volumes
MeshSize{ PointsOf{ Volume{:}; } } = lcar1;

// Override this constraint on the points of the cylinder
MeshSize{ PointsOf{ Volume{2}; } } = lcar3;

BooleanDifference{ Volume{1}; Delete; }{ Volume{2}; Delete; }

Coherence;

// Define the physical groups
Physical Surface(0) = {5};
Physical Surface(1) = {3};
Physical Surface(2) = {1,2,4,6};
Physical Surface(3) = {7,8,9,10,11,12};
Physical Volume(100) = {1};

// Generate the mesh
Mesh 3;

// Save the mesh
Save "cylinder3d-2.msh";
Save "cylinder3d-2.vtk";
