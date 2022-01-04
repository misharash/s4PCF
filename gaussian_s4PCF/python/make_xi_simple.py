import numpy as np
import os

r_low, r_high = np.loadtxt("radial_binning_corr.csv").T

norm = 50 
xi = 3*norm*(r_high - r_low)/(r_high**3 - r_low**3) # norm/r^2 averaged over bins
xi = xi[:, None] # add second dimension

r = (r_low + r_high)/2
mu = 0.5

dirname = "xi_simple"
os.makedirs(dirname, exist_ok=1)
filename = os.path.join(dirname, "xi_n200_m1_11")

header = " ".join([f"{rr:.8e}" for rr in r]) + "\n" + f"{mu:.8e}"
np.savetxt(filename, xi, fmt='%.8e', comments='', header=header)