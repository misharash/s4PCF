import os
import numpy as np

N = 50

header = None
arrays = []
for i in range(1, N+1):
    filename = f"xi_x50/{i:04d}/xi_n200_m120_11.dat"
    arrays.append(np.loadtxt(filename, skiprows=2))
    with open(filename) as f:
        header = ''.join(f.readlines()[:2]).strip()
arrays = np.array(arrays)
out = np.mean(arrays, axis=0)
dirname = f"xi_x{N}"
os.makedirs(dirname, exist_ok=1)
filename = os.path.join(dirname, "xi_n200_m120_11.dat")
np.savetxt(filename, out, fmt="%.8e", comments='', header=header)
