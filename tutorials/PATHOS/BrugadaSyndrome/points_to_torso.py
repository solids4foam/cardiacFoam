from paraview.simple import *

data = FindSource('ECG.foam')
view = GetActiveViewOrCreate('RenderView')


coords = [
    ("V1", [-26.1184,  -283.543,   -70.0558]),
    ("V2", [ 73.5820,  -278.089,   -70.1427]),
    ("V3", [ 85.1153,  -276.414,   -93.8405]),
    ("V4", [ 99.7976,  -269.695,  -106.5960]),
    ("V5", [122.8850,  -253.215,  -119.1550]),
    ("V6", [148.6130,  -227.840,  -119.3170]),
]


for name, pt in coords:
    probe = ProbeLocation(
        registrationName=name,
        Input=data,
        ProbeType='Fixed Radius Point Source'
    )
    probe.ProbeType.Center = pt

    display = Show(probe, view)
    display.SetRepresentationType('Points')
    display.PointSize = 50
    display.RenderPointsAsSpheres = 1

Render()
ResetCamera()
Render()