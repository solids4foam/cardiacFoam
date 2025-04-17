// Gmsh .geo file to create a truncated ellipsoid with an ellipsoid cavity
// corresponding to the idealised ventricle geometry in Land et al. (2015)

// Mesh spacing parameters
Include "meshSpacing.geo";

SetFactory("OpenCASCADE");

// Create outer sphere
Sphere(1) = {0,0,0,0.01};

// Create inner sphere
Sphere(2) = {0,0,0,0.007};

// Create box, which we will use for cutting out one half of the geometry
Box(3) = {-0.03, -0.03, -0.03, 0.06, 0.06, 0.035};

// Scale the outer sphere to be an ellipsoid
Dilate {{0,0,0}, {1,1,2}} { Volume{1}; }

// Scale the inner sphere to be an ellipsoid
Dilate {{0,0,0}, {1,1,2.4285714286}} { Volume{2}; }

// Subtract the inner ellipsoid from the outer ellipsoid
BooleanDifference{ Volume{1}; Delete; }{ Volume{2}; Delete; }

// Extract one eigth of the geometry
BooleanIntersection{ Volume{1}; Delete; }{ Volume{3}; Delete; }

// Tag the physical volume
Physical Volume("ventricle") = {1};

// Extract surfaces for boundary patches
// Get all boundary surfaces of the remaining volume
out[] = Boundary{ Volume{1}; };

// out[] now contains surfaces of:
// - inner cavity (inside)
// - outer wall (outside)
// - cut surface from the box (fixed)

// Use bounding boxes or IDs to distinguish and name the patches.
// This depends on your geometry, but assuming Gmsh assigns surface numbers in order,
// you can check in the GUI or mesh file which corresponds to what.
Physical Surface("inside") = {out[0]};
Physical Surface("outside") = {out[1]};
Physical Surface("fixed") = {out[2]};
