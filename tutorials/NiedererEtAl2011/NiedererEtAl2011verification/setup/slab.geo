// Niederer et al. (2011) slab as an unstructured tetrahedral mesh.
//
// Geometry is identical to the structured-hex verification case
// (../NiedererEtAl2011verification/system/blockMeshDict): a
// 20 x 3 x 7 mm slab, fibres along the long (x) axis. The only difference
// between this case and the hex one is HOW the mesh is built -- here gmsh
// generates near-uniform isotropic tets at characteristic length lc,
// converted with gmshToFoam. Everything downstream (conductivity, stimulus,
// ionic model, schemes, probes, reference activation times) is reused
// unchanged.
//
// The characteristic length placeholder 0.0001 is substituted by
// run_tet_sweep.sh (lc = the nominal dx in metres) to build the refinement
// ladder {0.5, 0.2, 0.1} mm, mirroring the published Niederer dx sweep.
//
// The box is built directly in METRES: gmshToFoam imports coordinates
// verbatim (no scale factor is applied), and cardiacFoam works in SI, so the
// hex case's blockMesh "scale 0.001" is folded into the dimensions here.
//
// A single physical surface ("walls") wraps all six faces. cardiacFoam
// creates Vm with a uniform zeroGradient (no-flux) boundary regardless of
// patch name -- the same insulated-slab boundary condition the hex case and
// the benchmark use -- so one patch is sufficient.

SetFactory("OpenCASCADE");

lc = 0.0001;

// 20 mm (x, fibre axis) x 3 mm (y) x 7 mm (z)
Box(1) = {0, 0, 0, 0.020, 0.003, 0.007};

Physical Volume("internal") = {1};
Physical Surface("walls") = {1, 2, 3, 4, 5, 6};

Mesh.MeshSizeMin = lc;
Mesh.MeshSizeMax = lc;
Mesh.Algorithm3D = 1;      // Delaunay: near-uniform, isotropic tets
Mesh.Optimize = 1;
Mesh.MshFileVersion = 2.2; // gmshToFoam requires the legacy msh2 format
