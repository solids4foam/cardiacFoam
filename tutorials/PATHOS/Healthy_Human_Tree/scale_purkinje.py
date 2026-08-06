import pyvista as pv

mesh = pv.read("purkinje.vtk")
mesh.points *= 0.001
mesh.save("purkinje_scaled.vtk", binary=False)
print("Scaled purkinje.vtk to purkinje_scaled.vtk (factor 0.001)")
