import sys
import os
from datetime import datetime
import numpy as np
from scipy.integrate import romberg
from scipy.interpolate import interp1d

if len(sys.argv) == 6 or len(sys.argv) == 7:
    xi_coarse_file = sys.argv[1]
    xi_fine_file = sys.argv[2]
    long_bin_file = sys.argv[3]
    short_bin_file = sys.argv[4]
    outfilename = sys.argv[5]
    if len(sys.argv) == 7:
        Nproc = int(sys.argv[6])
    else:
        Nproc = 1
else:
    raise Exception("Need to specify xi coarse file, xi fine file, long binning file, short binning file, output file name and (optional) number of processes!")

# read coarse xi in short bins
xi_bins, xi_bins_vals = np.loadtxt(xi_coarse_file)
xi_bins = xi_bins.astype(int)

# read fine xi data and create function
rmin = 1e-4
if os.path.exists(xi_fine_file):
    r_vals, xi_vals = np.loadtxt(xi_fine_file)
    r2xi_fun = interp1d(r_vals, r_vals**2 * xi_vals, kind='cubic', fill_value="extrapolate")
    xi_fun = lambda r: r2xi_fun(r) * np.fmax(r, rmin)**-2
else:
    r_vals = np.zeros(0)
    xi_fun = lambda r: np.fmax(r, rmin)**-2

# read bins
long_bins = np.loadtxt(long_bin_file)
short_bins = np.loadtxt(short_bin_file)
if len(xi_bins) != len(short_bins) or not all(xi_bins == np.arange(len(short_bins))):
    raise Exception("Short bins in xi coarse file don't match short bin file")

# dominant part - outer product of short bin-averaged CFs, independent of long bin
gs4PCF_base = np.ones(len(long_bins))[:, None, None] * np.expand_dims(np.outer(xi_bins_vals, xi_bins_vals), axis=0)

np.seterr(all='raise')

rel_precision = 1e-3 # order of desired relative precision for integrals
abs_precision = 1e-5 # order of desired absolute precision for integrals, should be fine as smallest complete integral is order of 1

def integration_wrapper(*args, **kwargs):
    return romberg(*args, rtol=rel_precision, tol=abs_precision, **kwargs)

def g(rij, R, rb, rc):
    F = lambda r: np.square(rb) * np.square(rc) * r - np.square(rb) * (r - R)**3 / 3 - np.square(rc) * (r - rij)**3 / 3 + np.square(rij) * np.square(R) * r - rij * R * (rij + R) * np.square(r) + (np.square(rij) + np.square(R) + 4 * rij * R) * r**3 / 3 - (R + rij) * r**4 / 2 + r**5 / 5
    rmin = np.maximum(rij - rb, R - rc)
    rmax = np.minimum(rij + rb, R + rc)
    return (F(rmax) - F(rmin)) * (rmax > rmin)

def G(rij, R, rb_min, rb_max, rc_min, rc_max):
    return g(rij, R, rb_max, rc_max) - g(rij, R, rb_min, rc_max) - g(rij, R, rb_max, rc_min) + g(rij, R, rb_min, rc_min)

def inner_integrand(rij, rb_min, rb_max, rc_min, rc_max):
    # xi_jk xi_il part - easier
    # xi_jk
    points_of_interest = np.array((rb_min, rb_max)) # throw r_bmin/max
    points_of_interest = rij + np.append(-points_of_interest, points_of_interest) # throw minus the abovementioned, add r_ij to all
    points_of_interest = np.append(points_of_interest, r_vals[np.logical_and(r_vals > min(points_of_interest), r_vals < max(points_of_interest))])
    # throw interpolation grid points fitting in the range
    points_of_interest = np.sort(np.unique(points_of_interest)) # delete repeating points and sort by ascension
    xi_jk_bavg = 0
    for (a, b) in zip(points_of_interest[:-1], points_of_interest[1:]):
        xi_jk_bavg += integration_wrapper(lambda R: R*xi_fun(R)*(np.square(rb_max) - np.fmax(np.square(rij - R), np.square(rb_min))), a, b, vec_func=True)
    xi_jk_bavg *= 3 / (4 * rij * (pow(rb_max, 3) - pow(rb_min, 3))) # common factor
    # xi_il
    points_of_interest = np.array((rc_min, rc_max)) # throw r_cmin/max
    points_of_interest = rij + np.append(-points_of_interest, points_of_interest) # throw minus the abovementioned, add r_ij to all
    points_of_interest = np.append(points_of_interest, r_vals[np.logical_and(r_vals > min(points_of_interest), r_vals < max(points_of_interest))])
    # throw interpolation grid points fitting in the range
    points_of_interest = np.sort(np.unique(points_of_interest)) # delete repeating points and sort by ascension
    xi_il_cavg = 0
    for (a, b) in zip(points_of_interest[:-1], points_of_interest[1:]):
        xi_il_cavg += integration_wrapper(lambda R: R*xi_fun(R)*(np.square(rc_max) - np.fmax(np.square(rij - R), np.square(rc_min))), a, b, vec_func=True)
    xi_il_cavg *= 3 / (4 * rij * (pow(rc_max, 3) - pow(rc_min, 3))) # common factor
    # carry their product to final result
    value = xi_jk_bavg * xi_il_cavg
    #print(f"Finished 1 part of inner integrand computation, rij={rij}, rb={rb_min, rb_max}, rc={rc_min, rc_max}") # comnment out later
    # xi_ij xi_kl part - harder
    points_of_interest = []
    for c1 in (rb_max, rb_min):
        for c2 in (rc_min, rc_max):
            for s1 in (-1, 1):
                for s2 in (-1, 1):
                    points_of_interest.append(s1 * c1 + s2 * c2) # throw points +/- r_bmin/max +/- r_cmin/max
    points_of_interest = rij + np.array(points_of_interest) # add r_ij to all above
    points_of_interest = np.append(points_of_interest, r_vals[np.logical_and(r_vals > min(points_of_interest), r_vals < max(points_of_interest))])
    # throw interpolation grid points fitting in the range
    points_of_interest = np.sort(np.unique(points_of_interest)) # delete repeating points and sort by ascension
    xi_ij = xi_fun(rij)
    xi_kl_bcavg = 0
    for (a, b) in zip(points_of_interest[:-1], points_of_interest[1:]):
        xi_kl_bcavg += integration_wrapper(lambda R: R*xi_fun(R)*G(rij, R, rb_min, rb_max, rc_min, rc_max), a, b, vec_func=True)
    xi_kl_bcavg *= 9 / (16 * rij * (pow(rb_max, 3) - pow(rb_min, 3)) * (pow(rc_max, 3) - pow(rc_min, 3))) # common factor
    value += xi_ij * xi_kl_bcavg
    return value * np.square(rij) # r_ij^2 weighting for bin average

def integrate_gs4PCF(ra_min, ra_max, rb_min, rb_max, rc_min, rc_max):
    points_of_interest = np.array((ra_min, ra_max))
    shifts = []
    for c1 in (rb_max, rb_min):
        for c2 in (rc_min, rc_max):
            for s1 in (-1, 0, 1):
                for s2 in (-1, 0, 1):
                    shifts.append(s1 * c1 + s2 * c2) # throw points +/-/0x r_bmin/max +/-/0x r_cmin/max, some loops are extra
    shifts = np.unique(shifts) # delete repetitions
    shited_r_vals = np.unique(np.ravel(r_vals[:, None] + shifts[None, :])) # these are potential discontinuities in derivatives, flattened and unique
    points_of_interest = np.sort(np.unique(np.append(points_of_interest, shited_r_vals[np.logical_and(shited_r_vals > ra_min, shited_r_vals < ra_max)]))) # throw array from above within integration limits, delete repetitions and sort
    value = 0
    for (a, b) in zip(points_of_interest[:-1], points_of_interest[1:]):
        value += integration_wrapper(inner_integrand, a, b, args=(rb_min, rb_max, rc_min, rc_max))
    return value * 3 / (pow(ra_max, 3) - pow(ra_min, 3)) # normalize by integral of r^2 in bin

def gs4PCF_integral_wrapper(i):
    ra_min, ra_max = long_bins[i]
    gs4PCF_i = np.zeros((len(short_bins), len(short_bins)))
    print(f"Started {i+1} of {len(long_bins)} ({datetime.now()})")
    for j, (rb_min, rb_max) in enumerate(short_bins):
        for k, (rc_min, rc_max) in enumerate(short_bins):
            if k < j: continue
            gs4PCF_i[j, k] = integrate_gs4PCF(ra_min, ra_max, rb_min, rb_max, rc_min, rc_max)
            gs4PCF_i[k, j] = gs4PCF_i[j, k] # symmetry
    print(f"Finished {i+1} of {len(long_bins)} ({datetime.now()})")
    return gs4PCF_i

if Nproc > 1:
    from multiprocessing import Pool
    pool = Pool(Nproc)
    print(f"Running with {Nproc} parallel processes")
    gs4PCF_integral = np.array(pool.map(gs4PCF_integral_wrapper, range(len(long_bins))))
elif Nproc == 1:
    print(f"Running sequentially")
    gs4PCF_integral = np.array([gs4PCF_integral_wrapper(i) for i in range(len(long_bins))])
else:
    raise Exception("Number of processes has to be positive!")
gs4PCF = gs4PCF_base + gs4PCF_integral

np.save(outfilename, gs4PCF)