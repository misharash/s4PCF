## Script to average xi from many data/mock files
## So far supports single field only

import sys
import os
import numpy as np

if len(sys.argv)!=5:
    print("Usage: python average_xi.py {WORKDIR} {N_R_BINS} {N_MU_BINS} {N_DATA}")
    sys.exit()
workdir = str(sys.argv[1])
nrbins = int(sys.argv[2])
nmubins = int(sys.argv[3])
ndata = int(sys.argv[4])

basename = f"xi_n{nrbins}_m{nmubins}_11.dat"

header = None
arrays = []
for i in range(1, ndata+1):
    filename = os.path.join(workdir, f"{i:04d}/{basename}")
    arrays.append(np.loadtxt(filename, skiprows=2)) # read data
    with open(filename) as f:
        new_header = ''.join(f.readlines()[:2]).strip() # read first two lines (have different format)
        if header is None:
            header = new_header # save first header
        elif new_header != header: # check that header does not change suddenly
            print(f"Warning: header (first two lines) of {filename} differs from the first header")
arrays = np.array(arrays)
out = np.mean(arrays, axis=0) # output simple average of corrfuncion, assumes about equal size of data/mock files
dirname = os.path.join(workdir, f"avg_{ndata}")
os.makedirs(dirname, exist_ok=1)
filename = os.path.join(dirname, basename)
np.savetxt(filename, out, fmt="%.8e", comments='', header=header) # save to file, don't append comment marks to header
