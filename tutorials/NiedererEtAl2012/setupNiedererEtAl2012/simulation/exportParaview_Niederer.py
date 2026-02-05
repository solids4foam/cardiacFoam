import sys
import os
from paraview.simple import *
import argparse
import os

def float_to_str(value):
    # Convert float to string without decimal dot, preserving precision as needed
    s = f"{value:.7f}"  # 6 decimals, e.g. 0.005000
    s = s.rstrip('0').rstrip('.')  # remove trailing zeros and dot
    s = s.replace('.', '')  # remove decimal point entirely
    return s
  

def extract_activation_data(case_path, points_file, output_dir, dx=None, dt=None,tissue = None, ionicModel = None, solver = None):
    #Define strings to pass for writing my csv files. 0.001 = 0001
    dt_str = float_to_str(dt)
    dx_str = float_to_str(dx)

    """
    Extracts activation data from an OpenFOAM case using ParaView (pvpython).
    
    Parameters:
    - case_path: Path to the .foam file (e.g. "/path/to/case.foam")
    - points_file: Path to Niederer points TXT file
    - output_dir: Directory where output CSVs will be saved
    - dx, dt: Optional metadata to include in output filenames
    """
    paraview.simple._DisableFirstRenderCameraReset()

    

    # Load OpenFOAM case
    casefoam = OpenFOAMReader(registrationName='casefoam', FileName=case_path)
    casefoam.MeshRegions = ['internalMesh']
    casefoam.CellArrays = ['activationTime']
    casefoam.UpdatePipeline()

    # Get final timestep
    timesteps = casefoam.TimestepValues
    last_time = timesteps[-1] if timesteps else 0.0

    # -------------------------------
    # Plot Over Line
    plotLine = PlotOverLine(Input=casefoam)
    plotLine.Point1 = [0.0, 0.0, 0.007]
    plotLine.Point2 = [0.02, 0.003, 0.0]
    plotLine.UpdatePipeline(time=last_time)

    plotLinePass = PassArrays(Input=plotLine)
    plotLinePass.PointDataArrays = ['activationTime', 'arc_length']
    plotLinePass.UpdatePipeline(time=last_time)

    # Output file for line plot
    line_output = os.path.join(output_dir, f"{solver}_{ionicModel}_{tissue}_line_DT{dt_str}_DX{dx_str}.csv" if dx and dt else "plot_line_output.csv")
    SaveData(line_output, proxy=plotLinePass)

    # -------------------------------
    # Read Niederer points
    pointsCSV = CSVReader(FileName=[points_file])
    pointsCSV.HaveHeaders = 1
    pointsCSV.FieldDelimiterCharacters = ' '
    pointsCSV.UpdatePipeline()

    points3D = TableToPoints(Input=pointsCSV)
    points3D.XColumn = 'X'
    points3D.YColumn = 'Y'
    points3D.ZColumn = 'Z'
    points3D.UpdatePipeline()

    # Resample
    resampled = ResampleWithDataset(SourceDataArrays=casefoam, DestinationMesh=points3D)
    resampled.UpdatePipeline(time=last_time)

    resampledPass = PassArrays(Input=resampled)
    resampledPass.PointDataArrays = ['Label', 'activationTime']
    resampledPass.UpdatePipeline(time=last_time)

    resample_output_points = os.path.join(output_dir, f"{solver}_{ionicModel}_{tissue}_points_DT{dt_str}_DX{dx_str}.csv" if dx and dt else "resample_output.csv")
    SaveData(resample_output_points, proxy=resampledPass)

    print(f"✅ Exported:\n - {line_output}\n - {resample_output_points}")


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Extract activation data with ParaView")
    parser.add_argument('--case', required=True, help='Path to .foam case file')
    parser.add_argument('--points', required=True, help='Path to points TXT file')
    parser.add_argument('--outdir', required=True, help='Output directory')
    parser.add_argument('--dx', type=float, required=False, help='DX value (mm)')
    parser.add_argument("--dt", type=float, required=False, help="Time step in ms")
    parser.add_argument('--tissue', default='test', help='Tissue model for filenames')
    parser.add_argument('--ionicModel', required=True, help='Ionic model name for filenames')
    parser.add_argument('--solver', required=False, help='solver type')

    args = parser.parse_args()

    extract_activation_data(
        case_path=args.case,
        points_file=args.points,
        output_dir=args.outdir,
        dx=args.dx,
        dt=args.dt,
        tissue= args.tissue,
        ionicModel=args.ionicModel,
        solver=args.solver
    )






