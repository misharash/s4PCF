### nbody_to_xyzw.py (Michael Rashkovetskyi, with the help of Boryana Hadzhiyska, 2022)
# This reads an ECSV file (as at least in AbacusSummit mocks) and conerts it to xyzw plain text format assumed by s4PCF code & running script
# Output contains real-space (not redshift-space) coordinates

import sys
import numpy as np
from astropy.io import ascii

if len(sys.argv) != 4:
    raise Exception("Need to specify input file, box size and output file!")
else:
    infile = sys.argv[1]
    boxsize = float(sys.argv[2])
    outfile = sys.argv[3]

gals_arr = ascii.read(infile) # contains x, y, z and so on, column titles specified in first non-comment string
# stack 3 coordinate columns and weights of 1 for each galaxy, transpose
gals_pos_weights = np.vstack((gals_arr['x'], gals_arr['y'], gals_arr['z'], np.ones_like(gals_arr['x'])))
gals_pos_weights[:3] += boxsize/2 # shift coordinates from [-boxsize/2, boxsize/2] range to [0, boxsize]
np.savetxt(outfile, gals_pos_weights.T) # save to file