import gmsh
import argparse
import math

parser = argparse.ArgumentParser()
parser.add_argument("--h", type=float, default=0.05)
parser.add_argument("--output", default="heart_iso.msh")
args = parser.parse_args()

gmsh.initialize()
gmsh.merge("meshes/biv_ellipsoid_surface.msh")

# Reconstruct geometry to allow surface remeshing
angle = math.pi / 4.0
gmsh.model.mesh.classifySurfaces(angle, True, True, angle)
gmsh.model.mesh.createGeometry()

# Create volume
surfaces = gmsh.model.getEntities(2)
sl = gmsh.model.geo.addSurfaceLoop([s[1] for s in surfaces])
gmsh.model.geo.addVolume([sl])
gmsh.model.geo.synchronize()

# Assign physical groups dynamically by bounding box
endo_rv_tags = []
for s in gmsh.model.getEntities(2):
    tag = s[1]
    bbox = gmsh.model.getBoundingBox(2, tag)
    # bbox is (xmin, ymin, zmin, xmax, ymax, zmax)
    xmax = bbox[3]
    if xmax < 1e-6:
        # BASE
        gmsh.model.addPhysicalGroup(2, [tag], 1)
        gmsh.model.setPhysicalName(2, 1, "BASE")
    elif xmax > 3.5:
        # EPI
        gmsh.model.addPhysicalGroup(2, [tag], 2)
        gmsh.model.setPhysicalName(2, 2, "EPI")
    elif xmax > 2.6 and xmax < 3.2:
        # ENDO_RV
        endo_rv_tags.append(tag)
    elif xmax > 2.0 and xmax <= 2.6:
        # ENDO_LV
        gmsh.model.addPhysicalGroup(2, [tag], 4)
        gmsh.model.setPhysicalName(2, 4, "ENDO_LV")

if endo_rv_tags:
    gmsh.model.addPhysicalGroup(2, endo_rv_tags, 3)
    gmsh.model.setPhysicalName(2, 3, "ENDO_RV")

# Assign physical group to volume
gmsh.model.addPhysicalGroup(3, [sl], 100)
gmsh.model.setPhysicalName(3, 100, "internal")


gmsh.option.setNumber("Mesh.CharacteristicLengthMax", args.h)
gmsh.option.setNumber("Mesh.CharacteristicLengthMin", args.h / 10.0)
# Output version 2 for gmshToFoam
gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)

gmsh.model.mesh.generate(3)
gmsh.write(args.output)
gmsh.finalize()
