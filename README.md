# s4PCF - code to compute the 4-point correllation function (4PCF) binned by distances between 2 different pairs of galaxies
Currently consists of two different executables:

* `s4PCF` in main directory - counts quadruplets (also pairs and triplets) in bins. Based on [encore](https://github.com/oliverphilcox/encore) by Oliver Philcox. `run_npcf.py` script (modernized version of `run_npcf.csh`) combines data and random catalogs, invokes the C++ code to count quadruplets, and divides counts to get 4,2,3-point correlation functions.
* `gs4PCF` in `gaussian_s4PCF` subdirectory - randomly throws quadruplets from random catalog and averages the disconnected/Gaussian contribution to 4PCF. Based on [RascalC](https://github.com/oliverphilcox/RascalC) by Oliver Philcox, comes with a bunch of useful python scripts.

Their better integration is work in progress.

## Binning
Quadruplets of galaxies A, B, C, D (or $i$, $j$, $k$, $l$) are binned by lengths of AB, AC, **BD** ($r_{ij}$, $r_{ik}$, $r_{jl}$) sides (isotropically), in contrast with AB, AC, **AD** ($r_{ij}$, $r_{ik}$, $r_{il}$) in `encore` (which also uses higher order multipoles).
The aim is to isolate two close pairs better.

AB ($r_{ij}$) is binned differently from two others, `NBIN_LONG` bins are uniformly spaced between `rmin_long` and `rmax_long`.

AC and BD ($r_{ik}$ and $r_{jl}$) are binned in the same way, `NBIN` bins are uniformly spaced between `rmin` and `rmax`.
BD ($r_{jl}$) bin number is more or equal than AC ($r_{ik}$) bin number, i.e. equal bins are included, and all bins can be recovered by symmetry.

`s4PCF` also outputs isotropic 2- and 3-point correlation functions in `NBIN` bins uniformly spaced between `rmin` and `rmax`, as a by-product.

## Self-counting and other concerns
If `rmin_long < 2 * rmax`, a triple loop is activated to exclude self-counting triangles (C=D or $k=l$), while everything else is done in double loops (much faster, but often still slow).

Additionally, when `rmin_long < rmax` other possibilities arise: A=D ($i=l$) or/and B=C ($j=k$), although these can be excluded using only a double loop.

Both codes now have an option (`PREVENT_TRIANGLES` in `s4PCF` and `prevent_triangles` in `gs4PCF`, set in `gaussian_s4PCF/modules/parameters.h`) to exclude bins that allow triangular configurations (AC+BD>AB, $r_{ik}+r_{jl}>r_{ij}$) and thus skip all the self-counting exclusion, although the results then may be tricky to read and interpret.

You might want to set `rmin_long >= 3 * rmax`, because then it is guaranteed that AC and BD ($r_{ik}$ and $r_{jl}$) are the shortest distances, and no self-counting issues will appear.
