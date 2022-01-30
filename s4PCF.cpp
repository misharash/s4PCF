// s4PCF.cpp -- Michael Rashkovetskyi, 2021. Based on encore by Oliver Philcox.

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/time.h>

#include "STimer.cc"
#include "threevector.hh"

// For multi-threading:
#ifdef OPENMP
#include <omp.h>
#endif

// NBIN is the number of bins we'll sort the radii into.
#define NBIN_SHORT 15 // short sides
#define NBIN_LONG 18 // long side
#define NBIN_CF 450 // fine, anisotropic 2-point correlation function
// MBIN is number of bins for mu
#define MBIN_CF 10 // fine, anisotropic 2-point correlation function

// Whether to exclude bins that can allow triangles (k=l), r_ij<=r_ik+r_jl
// Beneficial for performance - avoids triple loop
// Also guarantees no other 4PCF self-counts are involved
#define PREVENT_TRIANGLES 1

// MAXTHREAD is the maximum number of allowed threads.
// Big trouble if actual number exceeds this!
// No problem if actual number is smaller.
#define MAXTHREAD 40

typedef unsigned long long int uint64;

// Could swap between single and double precision here.
// Only double precision has been tested.
// Note that the AVX multipole code is always double precision.
typedef double Float;
// typedef float Float;
typedef double3 Float3;

// We need a vector floor3 function
Float3 floor3(float3 p) {
    return Float3(floor(p.x), floor(p.y), floor(p.z));
}

// we need a vector ceil3 function
Float3 ceil3(float3 p) {
    return Float3(ceil(p.x), ceil(p.y), ceil(p.z));
}

#define PAGE 4096  // To force some memory alignment.

// Classes specifying cells and grids
#include "modules/Basics.h"

// Pair counts class
#include "modules/Pairs.h"

Pairs pairs[MAXTHREAD];

#include "modules/FinePairs.h"

FinePairs finepairs[MAXTHREAD];

// Here's a simple structure for our normalized differences of the positions
typedef struct Xdiff {
    Float dx, dy, dz, w;
} Xdiff;

// Include the NPCF class here
#include "modules/NPCF.h"

NPCF npcf[MAXTHREAD];

void set_npcf(Float rmin_short, Float rmax_short, Float rmin_long, Float rmax_long) {
    for (int t = 0; t < MAXTHREAD; t++) {
        npcf[t].reset();
        npcf[t].calc_4pcf_indices(rmin_short, rmax_short, rmin_long, rmax_long);
    }
}

void sum_power() {
    // Just add up all of the threaded power into the zeroth element
    for (int t = 1; t < MAXTHREAD; t++) {
        npcf[0].sum_power(npcf + t);
        pairs[0].sum_power(pairs + t);
        finepairs[0].sum_power(finepairs + t);
    }
    return;
}

// Include class which creates the pairs
#include "modules/ComputePairs.h"

// Include class which creates / reads in particles and assigns them to a grid
#include "modules/Driver.h"

// ================================ main() =============================

void usage() {
    fprintf(stderr, "\nUsage for s4PCF:\n");
    fprintf(stderr,
            "   -in <file>: The input file (space-separated x,y,z,w).  "
            "Default sample.dat.\n");
    fprintf(stderr,
            "   -outstr <outstring>: String to prepend to the output "
            "file.  Default sample.\n");
    fprintf(stderr,
            "   -def: This allows one to accept the defaults without "
            "giving other entries.\n");
    fprintf(stderr,
            "   -rmin_short <rmin_short>: The minimum radius of the short "
            "side bin.  Default 0.\n");
    fprintf(stderr,
            "   -rmax_short <rmax_short>: The maximum radius of the short "
            "side bin.  Default 30.\n");
    fprintf(stderr,
            "   -rmin_long <rmin_long>: The minimum radius of the long "
            "side bin.  Default 30.\n");
    fprintf(stderr,
            "   -rmax_long <rmax_long>: The maximum radius of the long "
            "side bin.  Default 120.\n");
    fprintf(stderr,
            "   -rmin_cf <rmin_cf>: The minimum radius of the fine, anisotropic 2PCF bin.  Default 0.\n");
    fprintf(stderr,
            "   -rmax_cf <rmax_cf>: The maximum radius of the fine, anisotropic 2PCF bin.  Default 180.\n");
    fprintf(stderr, "\n");
    fprintf(stderr,
            "   -ran <np>: Ignore any file and use np random perioidic "
            "points instead.\n");
    fprintf(stderr,
            "   -box <boxsize> : The periodic size of the computational "
            "domain, if particles are thrown randomly.  Default 400.\n");
    fprintf(stderr,
            "   -scale <rescale>: How much to dilate the input "
            "positions by.  Default 1.\n");
    fprintf(stderr,
            "             Negative values causes =boxsize, rescaling "
            "unit cube to full periodicity\n");
    fprintf(stderr,
            "   -nside <nside>: The grid size for accelerating the "
            "pair count.  Default 8.\n");
    fprintf(stderr,
            "             Recommend having several grid cells per rmax_short.\n");
    fprintf(stderr, "\n");
    fprintf(stderr,
            "Other important parameters can only be set during "
            "compilations:\n");
    fprintf(stderr, "   NBIN_SHORT:  The number of radial bins for short sides.\n");
    fprintf(stderr, "   NBIN_LONG:  The number of radial bins for long side.\n");
    fprintf(stderr, "   NBIN_CF:  The number of radial bins for fine, anisotropic 2PCF.\n");
    fprintf(stderr, "   MBIN_CF:  The number of angular bins for fine, anisotropic 2PCF.\n");
    fprintf(stderr, "Similarly, the radial and mu bin spacings (currently linear) are hard-coded.\n");
    fprintf(stderr, "\n");
    fprintf(stderr,
            "    -balance: Rescale the negative weights so that the "
            "total weight is zero.\n");
    fprintf(stderr, "    -invert: Multiply all the weights by -1.\n");

    exit(1);
    return;
}

int main(int argc, char* argv[]) {
    // Important variables to set!  Here are the defaults:
    Float boxsize = 400;
    // The periodicity of the position-space cube. (overwritten if reading from
    // file)
    Float rescale = 1.0;  // If left zero or negative, set rescale=boxsize
    // The particles will be read from the unit cube, but then scaled by
    // boxsize.
    Float rmax_short = 30;
    // The maximum radius of the largest bin.
    Float rmin_short = 0;
    // The minimum radius of the smallest bin.
    Float rmax_long = 120;
    // The maximum radius of the long side bin.
    Float rmin_long = 30;
    // The minimum radius of the long side bin.
    Float rmax_cf = 180;
    // The maximum radius of fine 2pcf bin.
    Float rmin_cf = 0;
    // The minimum radius of fine 2pcf bin.
    int nside = 50;
    // The grid size, which should be tuned to match boxsize and rmax_short.
    // Don't forget to adjust this if changing boxsize!
    int make_random = 0;
    // If set, we'll just throw random periodic points instead of reading the
    // file
    int np = -1;  // Will be number of particles in a random distribution,
                  // but gets overwritten if reading from a file.
    int qbalance = 0, qinvert = 0;
    const char default_fname[] = "sample.dat";
    const char default_outstr[] = "sample";
    char* fname = NULL;
    char* outstr = NULL;

    // The periodicity of the position-space cuboid in 3D.
    Float3 rect_boxsize = {boxsize, boxsize,
                           boxsize};  // this is overwritten on particle read-in
    Float cellsize;

    STimer TotalTime, Prologue, Epilogue, PairTime, IOTime;
    // Detailed timings
    STimer InfileReadTime, WeightsReadTime, GridTime, OutputTime;

    TotalTime.Start();
    Prologue.Start();
    if (argc == 1)
        usage();
    int i = 1;
    while (i < argc) {
        if (!strcmp(argv[i], "-boxsize") || !strcmp(argv[i], "-box")) {
            Float tmp_box = atof(argv[++i]);
            rect_boxsize = {tmp_box, tmp_box, tmp_box};
        } else if (!strcmp(argv[i], "-rescale") || !strcmp(argv[i], "-scale"))
            rescale = atof(argv[++i]);
        else if (!strcmp(argv[i], "-rmax_short") || !strcmp(argv[i], "-max_short"))
            rmax_short = atof(argv[++i]);
        else if (!strcmp(argv[i], "-rmin_short") || !strcmp(argv[i], "-min_short"))
            rmin_short = atof(argv[++i]);
        else if (!strcmp(argv[i], "-rmax_long") || !strcmp(argv[i], "-max_long"))
            rmax_long = atof(argv[++i]);
        else if (!strcmp(argv[i], "-rmin_long") || !strcmp(argv[i], "-min_long"))
            rmin_long = atof(argv[++i]);
        else if (!strcmp(argv[i], "-rmax_cf") || !strcmp(argv[i], "-max_cf"))
            rmax_cf = atof(argv[++i]);
        else if (!strcmp(argv[i], "-rmin_cf") || !strcmp(argv[i], "-min_cf"))
            rmin_cf = atof(argv[++i]);
        else if (!strcmp(argv[i], "-nside") || !strcmp(argv[i], "-ngrid") ||
                 !strcmp(argv[i], "-grid"))
            nside = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-in"))
            fname = argv[++i];
        else if (!strcmp(argv[i], "-outstr"))
            outstr = argv[++i];
        else if (!strcmp(argv[i], "-balance"))
            qbalance = 1;
        else if (!strcmp(argv[i], "-invert"))
            qinvert = 1;
        else if (!strcmp(argv[i], "-ran") || !strcmp(argv[i], "-np")) {
            double tmp;
            if (sscanf(argv[++i], "%lf", &tmp) != 1) {
                fprintf(stderr, "Failed to read number in %s %s\n", argv[i - 1],
                        argv[i]);
                usage();
            }
            np = tmp;
            make_random = 1;
        } else if (!strcmp(argv[i], "-def") || !strcmp(argv[i], "-default")) {
            fname = NULL;
        }
        else {
            fprintf(stderr, "Don't recognize %s\n", argv[i]);
            usage();
        }
        i++;
    }

    // Compute smallest and largest boxsizes
    Float box_min = fmin(fmin(rect_boxsize.x, rect_boxsize.y), rect_boxsize.z);
    Float box_max = fmax(fmax(rect_boxsize.x, rect_boxsize.y), rect_boxsize.z);

    assert(i == argc);  // For example, we might have omitted the last
                        // argument, causing disaster.

    assert(box_min > 0.0);
    assert(rmax_short > 0.0);
    assert(rmin_short >= 0.0);
    assert(rmax_long > 0.0);
    assert(rmin_long >= 0.0);
    assert(nside > 0);
    assert(nside < 300);  // Legal, but rather unlikely that we should use
                          // something this big!
    if (rescale < 0.0)
        rescale = box_max;  // This would allow a unit cube to fill the periodic
                            // volume
    if (rescale == 0.0)
        rescale = 1;  // no rescaling
    if (fname == NULL)
        fname = (char*)default_fname;  // No name was given
    if (outstr == NULL)
        outstr = (char*)default_outstr;  // No outstring was given

    // Output for posterity
    printf("\nBox Size = {%6.5e,%6.5e,%6.5e}\n", rect_boxsize.x, rect_boxsize.y,
           rect_boxsize.z);
    printf("Grid = %d\n", nside);
    printf("Minimum Radius for short sides = %6.3g\n", rmin_short);
    printf("Maximum Radius for short sides = %6.3g\n", rmax_short);
    printf("Minimum Radius for long side = %6.3g\n", rmin_long);
    printf("Maximum Radius for long side = %6.3g\n", rmax_long);
    printf("Minimum Radius for fine 2-point correlation function = %6.3g\n", rmin_cf);
    printf("Maximum Radius for fine 2-point correlation function = %6.3g\n", rmax_cf);
    Float gridsize = rmax_short / (box_max / nside);
    printf("Max short radius in Grid Units = %6.3g\n", gridsize);
    if (gridsize < 1)
        printf("#\n# WARNING: grid appears inefficiently coarse\n#\n");
    printf("Bins = %d\n", NBIN_SHORT);

    IOTime.Start();
    InfileReadTime.Start();
    Particle* orig_p;
    Float3 shift;
    if (make_random) {
        // If you want to just make random particles instead:
        assert(np > 0);
        orig_p = make_particles(rect_boxsize, np);
        cellsize = rect_boxsize.x / nside;  // define size of cells
    } else {
        orig_p = read_particles(rescale, &np, fname);
        assert(np > 0);
        // Update boxsize here
        compute_bounding_box(orig_p, np, rect_boxsize, cellsize, fmax(rmax_short, rmax_long), shift, nside);
    }

    if (qinvert)
        invert_weights(orig_p, np);
    if (qbalance)
        balance_weights(orig_p, np);

    InfileReadTime.Stop();

    GridTime.Start();

    // Now ready to compute!
    // Sort the particles into the grid.
    Grid grid(orig_p, np, rect_boxsize, cellsize, shift);
    printf("# Done gridding the particles\n");
    printf("# %d particles in use, %d with positive weight\n", grid.np,
           grid.np_pos);
    printf("# Weights: Positive particles sum to %f\n", grid.sumw_pos);
    printf("#          Negative particles sum to %f\n", grid.sumw_neg);
    free(orig_p);

    Float grid_density = (double)np / grid.nf;
    printf("Average number of particles per grid cell = %6.2g\n", grid_density);
    printf("Average number of particles within allowed radii shell = %6.2g\n",
           np * 4.0 * M_PI / 3.0 * (pow(rmax_short, 3.0) - pow(rmin_short, 3.0)) /
               (rect_boxsize.x * rect_boxsize.y * rect_boxsize.z));
    if (grid_density < 1)
        printf("#\n# WARNING: grid appears inefficiently fine.\n#\n");

    GridTime.Stop();
    IOTime.Stop();

    set_npcf(rmin_short, rmax_short, rmin_long, rmax_long);
    fflush(NULL);

    Prologue.Stop();

    // Everything above here takes negligible time.  This line is nearly all of
    // the work.
    PairTime.Start();
    compute_pairs(&grid, rmin_short, rmax_short, rmin_long, rmax_long, rmin_cf, rmax_cf, np);
    printf("# Done counting the pairs\n");
    PairTime.Stop();

    // Output the results
    Epilogue.Start();
    OutputTime.Start();
    sum_power();
    OutputTime.Stop();

    // printf("\n# Binned weighted pair counts, monopole and quadrupole\n");
    // pairs[0].report_pairs();

    // Save the outputs
    pairs[0].save_pairs(outstr, rmin_short, rmax_short, grid.sumw_pos);
    finepairs[0].save_pairs(outstr, rmin_cf, rmax_cf, grid.sumw_pos);
    npcf[0].save_power(outstr, rmin_short, rmax_short, rmin_long, rmax_long, grid.sumw_pos);

    npcf[0].report_timings();

    Epilogue.Stop();
    TotalTime.Stop();
    printf("\n# Total Time: %4.1f s\n", TotalTime.Elapsed());
    printf("# Prologue: %6.3f s (%4.1f%%)\n", Prologue.Elapsed(),
           Prologue.Elapsed() / TotalTime.Elapsed() * 100.0);
    printf("# Epilogue: %6.3f s (%4.1f%%)\n", Epilogue.Elapsed(),
           Epilogue.Elapsed() / TotalTime.Elapsed() * 100.0);
    printf("# IO Time:  %6.3f s (%4.1f%%)\n", IOTime.Elapsed(),
           IOTime.Elapsed() / TotalTime.Elapsed() * 100.0);
    printf("# Pairs:    %6.3f s (%4.1f%%)\n", PairTime.Elapsed(),
           PairTime.Elapsed() / TotalTime.Elapsed() * 100.0);

    // Detailed timing breakdown
    printf("\n# Load Particles: %6.3f s (%4.1f%%)\n", InfileReadTime.Elapsed(),
           InfileReadTime.Elapsed() / TotalTime.Elapsed() * 100.0);
    printf("# Load Weights: %6.3f s (%4.1f%%)\n", WeightsReadTime.Elapsed(),
           WeightsReadTime.Elapsed() / TotalTime.Elapsed() * 100.0);
    printf("# Grid Allocation:  %6.3f s (%4.1f%%)\n", GridTime.Elapsed(),
           GridTime.Elapsed() / TotalTime.Elapsed() * 100.0);
    printf("# NPCF Output:    %6.3f s (%4.1f%%)\n", OutputTime.Elapsed(),
           OutputTime.Elapsed() / TotalTime.Elapsed() * 100.0);

    return 0;
}
