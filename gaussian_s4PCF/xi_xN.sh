#!/bin/bash

# define parameters
NTHREADS=10 # number of threads
NRBINS=900 # number of radial bins
NMUBINS=10 # number of angular (mu) bins
NQPM=50 # number of mock catalogues to use
OUTDIR=xi # output directory
MUMAX=1. # maximal mu
# cosmology
OMEGA_M=0.29 # Omega_m
OMEGA_K=0. # Omega_K
W_WCDM=-1 # w, here of cosmological constant

# download and uncompress SDSS DR12 randoms
wget https://data.sdss.org/sas/dr12/boss/lss/qpm_mocks/mock_random_DR12_CMASS_N_50x1.rdzw.gz
gunzip https://data.sdss.org/sas/dr12/boss/lss/qpm_mocks/mock_random_DR12_CMASS_N_50x1.rdzw.gz
# convert to xyz in Mpc/h
python python/convert_to_xyz.py mock_random_DR12_CMASS_N_50x1.rdzw qpm_randoms_50x.xyz ${OMEGA_M} ${OMEGA_K} ${W_WCDM}
# select small part randomly
python python/take_subset_of_particles.py qpm_randoms_50x.xyz qpm_randoms_10x.xyz 6420510
# 6420510 = 10 x 642051 (Ngalaxies in first mock)

# download SDSS DR12 QPM mocks archive, some 12 GB
wget https://data.sdss.org/sas/dr12/boss/lss/qpm_mocks/mock_galaxy_DR12_CMASS_N_QPM_allmocks.tar.gz
# create dir and unpack
mkdir -p mock_galaxy_DR12_CMASS_N_QPM
tar -xf mock_galaxy_DR12_CMASS_N_QPM_allmocks.tar.gz -C mock_galaxy_DR12_CMASS_N_QPM
# rename for convenience
for i in $(seq -f "%04.0f" ${NQPM}) ;
do
mv mock_galaxy_DR12_CMASS_N_QPM/{mock_galaxy_DR12_CMASS_N_QPM_,}${i}
done

mkdir -p qpm_galaxies

# convert rdz to xyz
for i in $(seq -f "%04.0f" ${NQPM}) ;
do
python python/convert_to_xyz.py mock_galaxy_DR12_CMASS_N_QPM/${i}.rdzw qpm_galaxies/${i}.xyz ${OMEGA_M} ${OMEGA_K} ${W_WCDM}
done

# do random counts once
PERIODIC=0 # aperiodic
NORMALIZE=1 # norm. by (sum w)^2
python python/RR_counts.py qpm_randoms_10x.xyzw radial_binning_corr.csv ${MUMAX} ${NMUBINS} ${NTHREADS} ${PERIODIC} ${OUTDIR}/ ${NORMALIZE}

# do correlation function
for i in $(seq -f "%04.0f" ${NQPM}) ;
do
python python/xi_estimator_aperiodic.py qpm_galaxies/${i}.xyz qpm_randoms_50x.xyz qpm_randoms_10x.xyz radial_binning_corr.csv ${MUMAX} ${NMUBINS} ${NTHREADS} ${OUTDIR}/${i}/ ${OUTDIR}/RR_counts_n${NRBINS}_m${NMUBINS}_11.txt;
done

# average xi
python python/average_xi.py ${OUTDIR} ${NRBINS} ${NMUBINS} ${NQPM}