# s4PCF - code to compute the 4-point correllation function binned by distances between 2 different pairs of galaxies
Based on [encore](https://github.com/oliverphilcox/encore) and [RascalC](https://github.com/oliverphilcox/RascalC) written by Oliver Philcox.

## Binning
Quadruplets of galaxies A, B, C, D (or $i$, $j$, $k$, $l$) are binned by lengths of AB, AC, **BD** ($r_{ij}$, $r_{ik}$, $r_{jl}$) sides (isotropically), in contrast with AB, AC, **AD** ($r_{ij}$, $r_{ik}$, $r_{il}$) in `encore` (which also uses higher order multipoles).
The aim is to isolate two close pairs better.

AB ($r_{ij}$) is binned differently from two others, `NBIN_LONG` bins are uniformly spaced between `rmin_long` and `rmax_long`.

AC and BD ($r_{ik}$ and $r_{jl}$) are binned in the same way, `NBIN` bins are uniformly spaced between `rmin` and `rmax`.
BD ($r_{jl}$) bin number is more or equal than AC ($r_{ik}$) bin number, i.e. equal bins are included, and all bins can be recovered by symmetry.

You might want to set `rmin_long >= 3 * rmax`, because then it is guaranteed that AC and BD ($r_{ik}$ and $r_{jl}$) are the shortest distances.
Also, if `rmin_long < 2 * rmax`, a triple loop is activated to exclude self-counting triangles (C=D or $k=l$), while everything else is done in double loops (much faster).
Additionally, when `rmin_long < rmax` other possibilities arise: A=D ($i=l$) or/and B=C ($j=k$), although these can be excluded using only a double loop.

The code also outputs isotropic 2- and 3-point correlation functions in `NBIN` bins uniformly spaced between `rmin` and `rmax`, as a by-product.
