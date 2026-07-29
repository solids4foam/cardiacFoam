import gmsh
import math

gmsh.initialize()
gmsh.open("/Users/simaocastro/noFrontendCardiacFoam_minor_errors/tutorials/FiberAlignedAnisotropicMesh/meshes/biv_ellipsoid_surface.msh")


gmsh.model.mesh.classifySurfaces(math.pi / 4.0, True, True, math.pi / 4.0)
gmsh.model.mesh.createGeometry()

surfaces = gmsh.model.getEntities(2)
for s in surfaces:
    tag = s[1]
    bbox = gmsh.model.getBoundingBox(2, tag)
    print(f"Surface {tag}: bbox {bbox}")
gmsh.finalize()
