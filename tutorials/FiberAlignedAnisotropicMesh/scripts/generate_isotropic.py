import gmsh
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--h", type=float, default=0.05)
parser.add_argument("--output", default="heart_iso.msh")
args = parser.parse_args()

gmsh.initialize()
gmsh.merge("meshes/biv_ellipsoid_surface.msh")

surfaces = [tag for dim, tag in gmsh.model.getEntities(2)]
sl = gmsh.model.geo.addSurfaceLoop(surfaces)
gmsh.model.geo.addVolume([sl])
gmsh.model.geo.synchronize()

gmsh.option.setNumber("Mesh.CharacteristicLengthMax", args.h)
gmsh.option.setNumber("Mesh.CharacteristicLengthMin", args.h / 10.0)
gmsh.option.setNumber("Mesh.Algorithm", 6)
gmsh.option.setNumber("Mesh.Algorithm3D", 10)

gmsh.model.mesh.generate(3)
gmsh.write(args.output)
gmsh.finalize()
