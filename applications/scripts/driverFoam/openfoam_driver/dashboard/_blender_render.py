"""Headless Blender script. Invoked as:
  blender --background --python _blender_render.py -- <surface> <out.blend> <out.png>

Imports a VTK/VTU/PLY/STL surface (converted to STL by render_case first if
needed), applies a simple red heart material, saves a .blend, and renders a PNG.
"""
import sys

import bpy


def _args():
    a = sys.argv
    return a[a.index("--") + 1:] if "--" in a else []


def main():
    surface, out_blend, out_png = _args()[:3]
    bpy.ops.wm.read_factory_settings(use_empty=True)
    if surface.lower().endswith(".stl"):
        try:
            bpy.ops.wm.stl_import(filepath=surface)
        except AttributeError:
            bpy.ops.import_mesh.stl(filepath=surface)
    elif surface.lower().endswith(".ply"):
        bpy.ops.wm.ply_import(filepath=surface)
    obj = bpy.context.selected_objects[0] if bpy.context.selected_objects else None

    mat = bpy.data.materials.new("heart")
    mat.use_nodes = True
    bsdf = mat.node_tree.nodes.get("Principled BSDF")
    if bsdf:
        bsdf.inputs["Base Color"].default_value = (0.72, 0.09, 0.11, 1.0)
    if obj is not None:
        obj.data.materials.append(mat)

    bpy.ops.object.light_add(type="SUN")
    bpy.ops.object.camera_add(location=(0, -3, 1.2), rotation=(1.2, 0, 0))
    bpy.context.scene.camera = bpy.context.object
    bpy.context.scene.render.filepath = out_png
    bpy.context.scene.render.image_settings.file_format = "PNG"
    bpy.ops.wm.save_as_mainfile(filepath=out_blend)
    bpy.ops.render.render(write_still=True)


if __name__ == "__main__":
    main()
