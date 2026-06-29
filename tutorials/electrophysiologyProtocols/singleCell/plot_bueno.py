import pandas as pd
import matplotlib.pyplot as plt
import os

base_dir = "/Users/simaocastro/noFrontendCardiacFoam_minor_errors/tutorials/coreProtocols/singleCell/postProcessing"

def plot_vm():
    plt.figure(figsize=(10, 6))
    
    files = [
        "BuenoOrovio_epicardialCells_S1_1000.txt",
        "BuenoOrovio_endocardialCells_S1_1000.txt",
        "BuenoOrovio_mCells_S1_1000.txt"
    ]
    
    for f in files:
        path = os.path.join(base_dir, f)
        if os.path.exists(path):
            df = pd.read_csv(path, sep=r'\s+')
            raw_name = f.replace("BuenoOrovio_", "").replace("_S1_1000.txt", "")
            plt.plot(df['time'], df['Vm'], label=raw_name)
            
    plt.xlabel("Time (s)")
    plt.ylabel("Vm (mV)")
    plt.title("Action Potential (Vm)")
    plt.legend()
    plt.grid(True)
    plt.savefig(os.path.join(base_dir, "plot_Vm.png"))
    plt.close()

def plot_ta():
    plt.figure(figsize=(10, 6))
    
    files = [
        "BuenoOrovio_epicardialCells_S1_1000_Ta.txt",
        "BuenoOrovio_endocardialCells_S1_1000_Ta.txt",
        "BuenoOrovio_mCells_S1_1000_Ta.txt"
    ]
    
    for f in files:
        path = os.path.join(base_dir, f)
        if os.path.exists(path):
            df = pd.read_csv(path, sep=r'\s+')
            raw_name = f.replace("BuenoOrovio_", "").replace("_S1_1000_Ta.txt", "")
            plt.plot(df['time'], df['Ta'], label=raw_name)
            
    plt.xlabel("Time (s)")
    plt.ylabel("Ta (kPa)")
    plt.title("Active Tension (Ta)")
    plt.legend()
    plt.grid(True)
    plt.savefig(os.path.join(base_dir, "plot_Ta.png"))
    plt.close()

if __name__ == "__main__":
    plot_vm()
    plot_ta()
