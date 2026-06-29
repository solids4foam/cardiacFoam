import pandas as pd
import numpy as np
import os

base_dir = "/Users/simaocastro/noFrontendCardiacFoam_minor_errors/tutorials/coreProtocols/singleCell/postProcessing"

def calc_apd90(file_name):
    path = os.path.join(base_dir, file_name)
    df = pd.read_csv(path, sep=r'\s+')
    
    # Only consider the first beat (time < 0.9s)
    df = df[df['time'] < 0.9]
    
    time = df['time'].values
    vm = df['Vm'].values
    
    # Calculate amplitude
    v_min = np.min(vm)
    v_max = np.max(vm)
    amplitude = v_max - v_min
    
    # 90% repolarization threshold
    v_repol90 = v_min + 0.1 * amplitude
    
    # Find upstroke time
    threshold = -10 # mV
    upstroke_idx = np.argmax(vm > threshold)
    t_upstroke = time[upstroke_idx]
    
    # Find repolarization time after peak
    peak_idx = np.argmax(vm)
    
    repol_idx = peak_idx + np.argmax(vm[peak_idx:] < v_repol90)
    t_repol = time[repol_idx]
    
    return t_repol - t_upstroke

try:
    apd_epi = calc_apd90("BuenoOrovio_epicardialCells_S1_1000.txt")
    print(f"APD90 'epicardialCells' (Bugged -> actually M-Cells): {apd_epi*1000:.2f} ms")
except Exception as e:
    print(f"Error epi: {e}")

try:
    apd_endo = calc_apd90("BuenoOrovio_endocardialCells_S1_1000.txt")
    print(f"APD90 'endocardialCells' (Correct -> Endo): {apd_endo*1000:.2f} ms")
except Exception as e:
    print(f"Error endo: {e}")

try:
    apd_m = calc_apd90("BuenoOrovio_mCells_S1_1000.txt")
    print(f"APD90 'mCells' (Bugged -> actually Epicardial): {apd_m*1000:.2f} ms")
except Exception as e:
    print(f"Error mCells: {e}")
