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
  

def extract_activation_data(case_path, output_dir, dx=None, dt=None,model = None, ionicModel = None):
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


    # Load the OpenFOAM case
    foam = OpenFOAMReader(registrationName='casefoam', FileName=case_path)
    foam.MeshRegions = ['internalMesh']
    foam.CellArrays = ['Vm', 'u1', 'u2', 'u3', 'ionicCurrent']
    foam.SkipZeroTime = 0
    foam.UpdatePipelineInformation()

    # Create a line sampling filter
    line = PlotOverLine(Input=foam)
    line.Point1 = [0.0, 0.05, 0.05]
    line.Point2 = [1.0, 0.05, 0.05]
    line.Resolution = 100  # number of sample points 
        # Output file for line plot
    # Output file for line plot
    line_output = os.path.join(
        output_dir, 
        f"{ionicModel}_{model}_{dx_str}_Cells.csv" if dx and dt else "plot_line_output.csv"
    )

    # Save all timesteps together in one CSV (includes Time column)
    SaveData(
        line_output,
        proxy=line,
        PointDataArrays=['Vm'],
        AddTime=1,        # ensures a "Time" column is written
        WriteTimeSteps=1  # ensures all timesteps are included
    )
    


    

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Extract activation data with ParaView")
    parser.add_argument('--case', required=True, help='Path to .foam case file')
    parser.add_argument('--outdir', required=True, help='Output directory')
    parser.add_argument('--dx', type=float, required=False, help='nCells')
    parser.add_argument("--dt", type=float, required=False, help="Time step in ms")
    parser.add_argument('--model', default='test', help='tissue Model name for filenames')
    parser.add_argument('--ionicModel', default='test', help='ionic Model name for filenames')

    args = parser.parse_args()

    extract_activation_data(
        case_path=args.case,
        output_dir=args.outdir,
        dx=args.dx,
        dt=args.dt,
        model=args.model,
        ionicModel= args.ionicModel,
    )





