### combine_files.py (Michael Rashkovetskyi, adapted from Oliver Philcox, 2021-2022)
# Handles both periodic and aperiodic cases
# In aperiodic case, this reads in a set of (data-random) and (random) particle counts and uses them to construct the N-point functions, including edge-correction
# In periodic case, this reads in only a set of (data-random) particle counts, computes analytic random counts and uses them to construct the N-point functions, without edge-correction
# It is designed to be used with the run_npcf.csh script
# Currently fine2PCF, 2PCF, 3PCF and 4PCF are supported.
# The output is saved to the working directory with the same format as the NPCF counts, with the filename ...zeta_{N}pcf.txt

import sys, os
import numpy as np

## First read-in the input file string from the command line
if len(sys.argv)!=3 or len(sys.argv)!=10:
    raise Exception("Need to specify the periodicity and input files!")
else:
    periodic = int(sys.argv[1])
    inputs = str(sys.argv[2])

if not periodic and len(sys.argv)!=3:
    raise Exception("In non-periodic case, need to specify only the periodicity (0) and input files!")

if periodic and len(sys.argv)!=10:
    raise Exception("In periodic case, need to specify the periodicity (1), input files, box size, minimal and maximal short radii, min and max long radii, min and max fine2PCF radii!")
else:
    boxsize = float(sys.argv[3])
    rmin_short = float(sys.argv[4])
    rmax_short = float(sys.argv[5])
    rmin_long = float(sys.argv[6])
    rmax_long = float(sys.argv[7])
    rmin_cf = float(sys.argv[8])
    rmax_cf = float(sys.argv[9])
    # number of galaxies is not needed, since counts are normalized by (sum of positive weights)^-N
    radius = lambda bin, maxbin, rmin, rmax: rmin + 1.*bin*(rmax-rmin)/(maxbin+1)
    bin_volume = lambda bin, maxbin, rmin, rmax: boxsize**-3*4.*np.pi/3.*(radius(bin+1,maxbin,rmin,rmax)**3.-radius(bin,maxbin,rmin,rmax)**3.)
    bin_volume_short = lambda bin, maxbin: bin_volume(bin, maxbin, rmin_short, rmax_short)
    bin_volume_long = lambda bin, maxbin: bin_volume(bin, maxbin, rmin_long, rmax_long)
    bin_volume_cf = lambda bin, maxbin: bin_volume(bin, maxbin, rmin_cf, rmax_cf)

print("Reading in files starting with %s\n"%inputs)

# Decide which N we're using
Ns = []
for N in ["2", "3", "4", "fine2"]:
    DmR_file = inroot+'.n00_%spcf.txt'%N
    if os.path.exists(DmR_file):
        Ns.append(N)

if len(Ns)==0:
    raise Exception("No files found with input string %s"%inputs)

for N in Ns:
    n = int(N[-1]) # same as N for 2,3,4; 2 for fine2

    if periodic:
        DmR_file = inroot+'.0.n0_%spcf.txt'%N
        bins = np.loadtxt(DmR_file, max_rows=n-1)
        if N == "fine2":
            bins = np.arange(len(bins))[:, None] # in this case "bins" is midpoints, need to override with bin numbers and add second dim
            n_mu = len(np.loadtxt(DmR_file, skiprows=1, max_rows=1)) # count mu bins
            countsR = 0.5 * bin_volume_cf(bins, np.max(bins)) / n_mu # factor of 1/2 because we count only i > j pairs for fine2pcf
        elif N == "2":
            countsR = bin_volume_short(bins, np.max(bins))
        elif N == "3":
            countsR = bin_volume_short(bins[0], np.max(bins)) * bin_volume_short(bins[1], np.max(bins))
        elif N == "4":
            countsR = bin_volume_long(bins[0], np.max(bins[0])) * bin_volume_short(bins[1], np.max(bins[1:])) * bin_volume_short(bins[2], np.max(bins[1:]))
        else: # should never reach this
            raise Exception(f"Unrecoginized N: {N}")
        R_file = inroot+'.n00_%spcf.txt'%N # technically it's DmR_file, but that will be overridden. R_file won't exist in periodic case. This will be used as source of bin data for output.
    else:
        # First load in R piece
        R_file = inputs+'.r_%spcf.txt'%N
        countsR = np.loadtxt(DmR_file, skiprows=(len(N)>1))[n-1:] # skipping rows with radial bins, skip 1 more row for fine2

    # Now load in D-R pieces and average
    countsN_all = []
    for i in range(100):
        DmR_file = inputs+'.n%s_%spcf.txt'%(str(i).zfill(2),N)
        if not os.path.exists(DmR_file): continue
        # Extract counts
        countsN_all.append(np.loadtxt(DmR_file, skiprows=(len(N)>1))[n-1:]) # skipping rows with radial bins, skip 1 more row for fine2
    countsN_all = np.asarray(countsN_all)
    N_files = len(countsN_all)
    countsN = np.mean(countsN_all,axis=0)
    # could use this next line to compute std from finite number of randoms!
    #countsNsig = np.std(countsN_all,axis=0)/np.sqrt(N_files)

    # Now compute edge-correction equations

    # isotropic 2/3/4PCF are all easy and similar
    zeta = countsN/countsR

    # Now save the output to file, copying the first few lines from the N files
    zeta_file = inputs+'.zeta_%spcf.txt'%N
    rfile = open(R_file,"r")
    zfile = open(zeta_file,"w")
    # copy comments and bins from first random file
    lnc = 0 # counter of lines that are not comments
    for line in rfile:
        if lnc >= n-1+(len(N)>1): break # only need N-1 data lines (2 for fine2PCF), i.e. radial bins, can terminate loop afterwards
        if line[0] != "#": lnc += 1
        zfile.write(line)
    rfile.close()
    # write NPCF
    for a in range(len(zeta)):
        for b in range(len(zeta[0])):
            zfile.write("%.8e\t" % zeta[a, b])
        zfile.write("\n")
    zfile.close()

    print("Computed %sPCF using %d (data-random) files, saving to %s\n"%(N,N_files,zeta_file))
