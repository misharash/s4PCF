## Script to average xi over mu bins
## So far supports single field only

import sys
import numpy as np

if len(sys.argv)!=3:
    print("Usage: python average_xi_radial.py {INFILE} {OUTFILE}")
    sys.exit()
infile = str(sys.argv[1])
outfile = str(sys.argv[2])

xi_in_data = np.loadtxt(infile, skiprows=2, ndmin=2)
r_bins = np.loadtxt(infile, max_rows=1)

xi_out_data = np.mean(xi_in_data, axis=1)
# RR counts in equally spaced mu bins are expected to be the same, thus simple averaging is correct

np.savetxt(outfile, np.array((r_bins, xi_out_data)), fmt="%.8e")
