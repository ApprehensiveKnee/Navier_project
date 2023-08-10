//+
SetFactory("OpenCASCADE");

lcar1 = 0.25;
lcar2 = .05;
lcar3 = .05;
// first generate a box with the given dimensions:
// 3,5 along x
// 2 along y
// 1 along z
// bottom left corner at (0,0,0)
Box(1) = {0,0,0, 5,3,1.5};

// then generate a cylinder with the given dimensions:
// radius 0.2
// height 1
// axis along y
// bottom center at (0.75,1,0.5)
Cylinder(2) = {0.75,1,0.75, 0,1,0, 0.2, 2*Pi};

// Assign a mesh size to all the points of all the volumes:
MeshSize{ PointsOf{ Volume{:}; } } = lcar1;

// Override this constraint on the points of the cylinder:
MeshSize{ PointsOf{ Volume{2}; } } = lcar3;


// Apply boolean operation to subtract the cylinder from the box
BooleanDifference(3) = { Volume{1};Delete; }{ Volume{2}; Delete; };

// Define physical surfaces for the box boundaries
Physical Volume(100) = {3};
Physical Surface(0) = {1};   // The surface on the positive x side of the box  --> Inlet
Physical Surface(1) = {6};    // The surface on the negative x side of the box --> Outlet
Physical Surface(2) = {2,3,4,5}; // The surface representing the hole --> Walls
// Define physical surface for the hole
Physical Surface(3) = {7,8,9}; // The surface representing the hole  ---> Cylinder

// Physical Surface("Right") = {2};   // The surface on the positive y side of the box
// Physical Surface("Left") = {4};    // The surface on the negative y side of the box
// Physical Surface("Bottom") = {5};  // The surface on the bottom of the box
// Physical Surface("Top") = {3};     // The surface on the top of the box


Coherence;

General.NumThreads = 0;
Geometry.OCCParallel = 1;
Geometry.OCCSewFaces = 1;

Mesh.Algorithm3D = 10;

// Finally, mesh the volume
Mesh 3;

// And save the mesh to file
Save "cylinder3d.msh";