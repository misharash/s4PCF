import sys
import os
from datetime import datetime
import numpy as np
from scipy.integrate import romberg
from scipy.interpolate import interp1d
from multiprocessing import Pool

if len(sys.argv) == 7:
    xi_coarse_file = sys.argv[1]
    xi_fine_file = sys.argv[2]
    long_bin_file = sys.argv[3]
    short_bin_file = sys.argv[4]
    outfilename = sys.argv[5]
    Nproc = int(sys.argv[6])
else:
    raise Exception("Need to specify xi coarse file, xi fine file, long binning file, short binning file, output file name and number of processes!")

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

def G(rij, R, rb_min, rb_max, rc_min, rc_max):
    F11 = lambda r: r * (np.square(rb_max) - np.square(rb_min)) * (np.square(rc_max) - np.square(rc_min))
    F12 = lambda r: (np.square(rb_max) - np.square(rb_min)) * (np.square(rc_max) * r - (r - R)**3 / 3)
    F21 = lambda r: (np.square(rc_max) - np.square(rc_min)) * (np.square(rb_max) * r - (r - rij)**3 / 3)
    F22 = lambda r: np.square(rb_max) * np.square(rc_max) * r - np.square(rb_max) * (r - R)**3 / 3 - np.square(rc_max) * (r - rij)**3 / 3 + np.square(rij) * np.square(R) * r - rij * R * (rij + R) * np.square(r) + (np.square(rij) + np.square(R) + 4 * rij * R) * r**3 / 3 - (R + rij) * r**4 / 2 + r**5 / 5
    limits1 = np.array((rij - rb_max, rij - rb_min, rij + rb_min, rij + rb_max))
    limits2 = np.array((R - rc_max, R - rc_min, R + rc_min, R + rc_max))
    integral_index = (1, 0, 1)
    integrals = ((F11, F12), (F21, F22))
    current = max(limits1[0], limits2[0])
    i = np.count_nonzero(limits1 <= current) - 1
    j = np.count_nonzero(limits2 <= current) - 1
    value = 0
    while i < len(limits1)-1 and j < len(limits2)-1:
        i_next, j_next = i, j
        if limits1[i + 1] <= limits2[j + 1]:
            i_next += 1
            nextel = limits1[i_next]
        else:
            j_next += 1
            nextel = limits2[j_next]
        integral = integrals[integral_index[i]][integral_index[j]]
        value += integral(nextel) - integral(current)
        i, j, current = i_next, j_next, nextel
    return value

def inner_integrand(rij, rb_min, rb_max, rc_min, rc_max):
    # xi_jk xi_il part - easier
    # xi_jk
    points_of_interest = np.array((rb_min, rb_max)) # throw r_bmin/max
    points_of_interest = rij + np.append(-points_of_interest, points_of_interest) # throw minus the abovementioned, add r_ij to all
    points_of_interest = np.append(points_of_interest, r_vals[np.logical_and(r_vals > min(points_of_interest), r_vals < max(points_of_interest))])
    # throw interpolation grid points fitting in the range
    points_of_interest = np.sort(np.unique(points_of_interest)) # delete repeating points and sort by ascension
    interval_len = points_of_interest[-1] - points_of_interest[0]
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
    interval_len = points_of_interest[-1] - points_of_interest[0]
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
    interval_len = points_of_interest[-1] - points_of_interest[0]
    xi_ik = xi_fun(rij)
    xi_kl_bcavg = 0
    for (a, b) in zip(points_of_interest[:-1], points_of_interest[1:]):
        xi_kl_bcavg += integration_wrapper(lambda R: R*xi_fun(R)*G(rij, R, rb_min, rb_max, rc_min, rc_max), a, b)
    xi_kl_bcavg *= 3 / (16 * rij * (pow(rb_max, 3) - pow(rb_min, 3)) * (pow(rc_max, 3) - pow(rc_min, 3))) # common factor
    value += xi_ik * xi_kl_bcavg
    return value * np.square(rij) # r_ij^2 weighting for bin average

def integrate_gs4PCF(ra_min, ra_max, rb_min, rb_max, rc_min, rc_max):
    interval_len = ra_max - ra_min
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

def integrate_gs4PCF_wrapper(i, j, k):
    print(f"Started {i + 1, j + 1, k + 1} of {len(long_bins), len(short_bins), len(short_bins)} ({datetime.now()})")
    value = integrate_gs4PCF(long_bins[i, 0], long_bins[i, 1], short_bins[j, 0], short_bins[j, 1], short_bins[k, 0], short_bins[k, 1])
    print(f"Finished {i + 1, j + 1, k + 1} of {len(long_bins), len(short_bins), len(short_bins)} ({datetime.now()})")
    return value

index_grid = np.meshgrid(np.arange(len(long_bins)), np.arange(len(short_bins)), np.arange(len(short_bins)))
index_grid = np.array([x.flatten() for x in index_grid])
index_grid = index_grid[:, index_grid[2] >= index_grid[1]] # k >= j in function

pool = Pool(Nproc)
gs4PCF_tmp_1D = np.array(pool.starmap(integrate_gs4PCF_wrapper, index_grid.T))
gs4PCF_integrated = np.zeros((len(long_bins), len(short_bins), len(short_bins)))
gs4PCF_integrated[index_grid[0], index_grid[1], index_grid[2]] = gs4PCF_tmp_1D
gs4PCF_integrated[index_grid[0], index_grid[2], index_grid[1]] = gs4PCF_tmp_1D # symmetry

gs4PCF = gs4PCF_base + gs4PCF_integrated

np.save(outfilename, gs4PCF)