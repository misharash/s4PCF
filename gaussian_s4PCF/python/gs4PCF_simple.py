import sys
import os
from datetime import datetime
import numpy as np
from scipy.integrate import romberg
from scipy.interpolate import interp1d

if len(sys.argv) == 6:
    xi_coarse_file = sys.argv[1]
    xi_fine_file = sys.argv[2]
    long_bin_file = sys.argv[3]
    short_bin_file = sys.argv[4]
    outfilename = sys.argv[5]
else:
    raise Exception("Need to specify xi coarse file, xi fine file, long binning file, short binning file and output file name!")

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
    xi_fun = lambda r: np.fmax(r, rmin)**-2

# read bins
long_bins = np.loadtxt(long_bin_file)
short_bins = np.loadtxt(short_bin_file)
if len(xi_bins) != len(short_bins) or not all(xi_bins == np.arange(len(short_bins))):
    raise Exception("Short bins in xi coarse file don't match short bin file")

# dominant part - outer product of short bin-averaged CFs, independent of long bin
gs4PCF = np.ones(len(long_bins))[:, None, None] * np.expand_dims(np.outer(xi_bins_vals, xi_bins_vals), axis=0)

np.seterr(all='raise')

precision = 1e-3 # order of desired precision for integrals

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
    xi_jk_bavg = romberg(lambda R: R*xi_fun(R)*(np.square(rb_max) - np.square(rij - R)), rij - rb_max, rij - rb_min, vec_func=True, rtol=precision*2*rb_max/(rb_max-rb_min))
    xi_jk_bavg += romberg(lambda R: R*xi_fun(R)*(np.square(rb_max) - np.square(rb_min)), rij - rb_min, rij + rb_min, vec_func=True, rtol=precision*rb_max/rb_min)
    xi_jk_bavg += romberg(lambda R: R*xi_fun(R)*(np.square(rb_max) - np.square(rij - R)), rij + rb_min, rij + rb_max, vec_func=True, rtol=precision*2*rb_max/(rb_max-rb_min))
    xi_jk_bavg *= 3 / (4 * rij * (pow(rb_max, 3) - pow(rb_min, 3))) # common factor
    # xi_il
    xi_il_cavg = romberg(lambda R: R*xi_fun(R)*(np.square(rc_max) - np.square(rij - R)), rij - rc_max, rij - rc_min, vec_func=True, rtol=precision*2*rc_max/(rc_max-rc_min))
    xi_il_cavg += romberg(lambda R: R*xi_fun(R)*(np.square(rc_max) - np.square(rc_min)), rij - rc_min, rij + rc_min, vec_func=True, rtol=precision*rc_max/rc_min)
    xi_il_cavg += romberg(lambda R: R*xi_fun(R)*(np.square(rc_max) - np.square(rij - R)), rij + rc_min, rij + rc_max, vec_func=True, rtol=precision*2*rc_max/(rc_max-rc_min))
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
                    points_of_interest.append(s1 * c1 + s2 * c2)
    points_of_interest = rij + np.sort(np.unique(points_of_interest)) # delete repeating points, sort and add rij
    interval_len = points_of_interest[-1] - points_of_interest[0]
    xi_ik = xi_fun(rij)
    xi_kl_bcavg = 0
    for (a, b) in zip(points_of_interest[:-1], points_of_interest[1:]):
        xi_kl_bcavg += romberg(lambda R: R*xi_fun(R)*G(rij, R, rb_min, rb_max, rc_min, rc_max), a, b, rtol=precision*interval_len/(b-a))
        # allow worse precision for short intervals
    xi_kl_bcavg *= 3 / (16 * rij * (pow(rb_max, 3) - pow(rb_min, 3)) * (pow(rc_max, 3) - pow(rc_min, 3))) # common factor
    value += xi_ik * xi_kl_bcavg
    return value * np.square(rij) # r_ij^2 weighting for bin average

def integrate_gs4PCF(ra_min, ra_max, rb_min, rb_max, rc_min, rc_max):
    value = romberg(inner_integrand, ra_min, ra_max, args=(rb_min, rb_max, rc_min, rc_max), rtol=precision)
    return value * 3 / (pow(ra_max, 3) - pow(ra_min, 3)) # normalize by integral of r^2 in bin

for i, (ra_min, ra_max) in enumerate(long_bins):
    print(f"Started {i+1} of {len(long_bins)} ({datetime.now()})")
    for j, (rb_min, rb_max) in enumerate(short_bins):
        for k, (rc_min, rc_max) in enumerate(short_bins):
            if k < j: continue
            print(f"Started {j + 1, k + 1} of {len(short_bins), len(short_bins)} ({datetime.now()})")
            gs4PCF[i, j, k] += integrate_gs4PCF(ra_min, ra_max, rb_min, rb_max, rc_min, rc_max)
            gs4PCF[i, k, j] = gs4PCF[i, j, k] # symmetry
            print(f"Finished {j + 1, k + 1} of {len(short_bins), len(short_bins)} ({datetime.now()})")
    print(f"Finished {i+1} of {len(long_bins)} ({datetime.now()})")

np.save(outfilename, gs4PCF)