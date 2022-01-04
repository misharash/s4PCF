import numpy as np
import os

r_low, r_high = np.loadtxt("radial_binning_corr.csv").T

nmu = 10

norm = 50 
xi = 3*norm*(r_high - r_low)/(r_high**3 - r_low**3) # norm/r^2 averaged over bins
xi = xi[:, None] # add second dimension
xi = np.repeat(xi, nmu, axis=1) # repeat for each mu bin

r = (r_low + r_high)/2
mu = np.linspace(0, 1, 2*nmu+1)[1::2]

dirname = "xi_simple"
os.makedirs(dirname, exist_ok=1)
filename = os.path.join(dirname, "xi_n200_m1_11.dat")

header = " ".join([f"{rr:.8e}" for rr in r]) + "\n" + " ".join([f"{mumu:.8e}" for mumu in mu])
np.savetxt(filename, xi, fmt='%.8e', comments='', header=header)