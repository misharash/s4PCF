### gen_periodic_randoms.py (Michael Rashkovetskyi, 2022)
# This generates uniform random positions within a cubic (periodic) box

import sys
import numpy as np

if len(sys.argv) != 4:
    raise Exception("Need to specify number of points, box size and output file!")
else:
    N = int(sys.argv[1])
    boxsize = float(sys.argv[2])
    outfile = sys.argv[3]

gals_pos = np.random.rand((3, N)) * boxsize # 3 coordinate sequences of length N, between 0 and boxsize
gals_pos_weights = np.append(gals_pos, np.ones((1, N)), axis=0) # append column of ones for weights
np.savetxt(outfile, gals_pos_weights.T) # save to file