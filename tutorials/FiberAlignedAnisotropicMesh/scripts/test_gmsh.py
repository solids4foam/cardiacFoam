import gmsh

gmsh.initialize()
gmsh.merge("meshes/biv_ellipsoid_surface.msh")
gmsh.merge("metric_test.pos")

bg_field = gmsh.model.mesh.field.add("PostView")
gmsh.model.mesh.field.setNumber(bg_field, "ViewIndex", 0)
gmsh.model.mesh.field.setAsBackgroundMesh(bg_field)

gmsh.option.setNumber("Mesh.Algorithm", 6)
gmsh.option.setNumber("Mesh.Algorithm3D", 1)
gmsh.option.setNumber("Mesh.AnisoMax", 10.0)
gmsh.option.setNumber("Mesh.CharacteristicLengthExtendFromBoundary", 0)

gmsh.model.mesh.generate(3)

print("Nodes:", gmsh.model.mesh.getNodes()[0].size)
gmsh.finalize()
