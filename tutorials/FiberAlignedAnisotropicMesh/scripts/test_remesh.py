import gmsh
import math

gmsh.initialize()
gmsh.merge("meshes/biv_ellipsoid_surface.msh")

# Extract edges and surfaces
angle = math.pi / 4.0
gmsh.model.mesh.classifySurfaces(angle, True, True, angle)
gmsh.model.mesh.createGeometry()

surfaces = gmsh.model.getEntities(2)
sl = gmsh.model.geo.addSurfaceLoop([s[1] for s in surfaces])
gmsh.model.geo.addVolume([sl])
gmsh.model.geo.synchronize()

gmsh.option.setNumber("Mesh.CharacteristicLengthMin", 0.04)
gmsh.option.setNumber("Mesh.CharacteristicLengthMax", 0.04)
gmsh.model.mesh.generate(3)

print("Surface nodes:", gmsh.model.mesh.getNodes(2, surfaces[0][1])[0].size)
gmsh.finalize()
