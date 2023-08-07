// Gmsh project created on Sun Jul 30 18:16:31 2023


// Define the characteristic length for mesh elements (adjust as needed)
lc = 0.1;

// Define points for the rectangular domain
Point(1) = {0, 0, 0, lc};
Point(2) = {3, 0, 0, lc};
Point(3) = {3, 1, 0, lc};
Point(4) = {0, 1, 0, lc};

// Define points for the circular hole (adjust the center and radius as needed)
Point(5) = {0.5, 0.5, 0, lc / 4}; // Center of the circle
Point(6) = {0.5 + 0.2, 0.5, 0, lc / 4}; // Point 1 on the circle
Point(7) = {0.5 , 0.5 + 0.2, 0, lc / 4}; // Point 1 on the circle
Point(8) = {0.5 - 0.2, 0.5, 0, lc / 4}; // Point 1 on the circle
Point(9) = {0.5, 0.5 - 0.2, 0, lc / 4}; // Point 1 on the circle

// Define lines for the rectangular domain
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

// Define a circle for the hole
Circle(5) = {6, 5, 7};
Circle(6) = {7, 5, 8};
Circle(7) = {8, 5, 9};
Circle(8) = {9, 5, 6};

// Define the loop for the outer boundary (rectangle)
Line Loop(10) = {1, 2, 3, 4};

// Define the loop for the circular hole
Line Loop(20) = {5, 6, 7, 8};

// Define the plane surface of the rectangular domain
Plane Surface(100) = {10, 20};

// Define the compound surface of the rectangular domain with a hole
//Compound Surface{100, 200};

// Define the physical surfaces (optional but useful for boundary conditions)
Physical Surface(100) = {100}; // Inlet = 1
Physical Curve(1) = {4}; // Inlet = 1
Physical Curve(2) = {2}; // Outlet = 2
Physical Curve(3) = {1}; // Bottom = 3
Physical Curve(4) = {3}; // Top = 4
Physical Curve(5) = {5,6,7,8}; // Cylinder = 5
Coherence;


// Mesh the compound surface
Mesh 2;

// Save the mesh
Save "cylinder2d.msh";
Save "cylinder2d.vtk";
