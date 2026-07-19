// Conformal tetrahedral bath--myocardium--bath mesh.
// 0.025000000000000001 is replaced by setup/run_mesh_gate.sh.

SetFactory("OpenCASCADE");

lc = 0.025000000000000001;
eps = 1e-7;

Box(1) = {-1, 0, 0, 1, 1, 1};
Box(2) = { 0, 0, 0, 1, 1, 1};
Box(3) = { 1, 0, 0, 1, 1, 1};

// Fragment all volumes together so the x=0 and x=1 interfaces share the
// same points and triangles. They are intentionally not Physical Surfaces:
// gmshToFoam must import them as internal faces, not boundary patches.
BooleanFragments{ Volume{1, 2, 3}; Delete; }{}

leftBath[] = Volume In BoundingBox{-1-eps, -eps, -eps, 0+eps, 1+eps, 1+eps};
heart[] = Volume In BoundingBox{0-eps, -eps, -eps, 1+eps, 1+eps, 1+eps};
rightBath[] = Volume In BoundingBox{1-eps, -eps, -eps, 2+eps, 1+eps, 1+eps};

Physical Volume("myocardium") = {heart[]};
Physical Volume("bath") = {leftBath[], rightBath[]};

xMin[] = Surface In BoundingBox{-1-eps, -eps, -eps, -1+eps, 1+eps, 1+eps};
xMax[] = Surface In BoundingBox{2-eps, -eps, -eps, 2+eps, 1+eps, 1+eps};
yMin[] = Surface In BoundingBox{-1-eps, -eps, -eps, 2+eps, eps, 1+eps};
yMax[] = Surface In BoundingBox{-1-eps, 1-eps, -eps, 2+eps, 1+eps, 1+eps};
zMin[] = Surface In BoundingBox{-1-eps, -eps, -eps, 2+eps, 1+eps, eps};
zMax[] = Surface In BoundingBox{-1-eps, -eps, 1-eps, 2+eps, 1+eps, 1+eps};

Physical Surface("xMin") = {xMin[]};
Physical Surface("xMax") = {xMax[]};
Physical Surface("sides") = {yMin[], yMax[], zMin[], zMax[]};

Mesh.MeshSizeMin = lc;
Mesh.MeshSizeMax = lc;
Mesh.Algorithm3D = 1;
Mesh.Optimize = 1;
Mesh.MshFileVersion = 2.2;
