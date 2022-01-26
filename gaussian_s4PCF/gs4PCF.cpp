// gs4PCF.cpp -- Michael Rashkovetskyi, 2021 based on grid_covariance.cpp of Oliver Philcox / Alexander Wiegand.

#include <sys/time.h>
#include <sys/resource.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <complex>
#include <algorithm>
#include <random>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>
#include "../threevector.hh"
#include <gsl/gsl_rng.h>
#include "ransampl/ransampl.h"
#include "../STimer.cc"
#include "cubature/cubature.h"
#include <limits>
#include <sys/stat.h>

// For multi-threading:
#ifdef OPENMP
#include <omp.h>
#include <sched.h>
#endif

// In order to not print the output matrices

#define PAGE 4096     // To force some memory alignment.

typedef unsigned long long int uint64;

// Could swap between single and double precision here.
typedef double Float;
typedef double3 Float3;


// Define module files
#include "modules/parameters.h"
#include "modules/cell_utilities.h"
#include "modules/grid.h"
#include "modules/correlation_function.h"
#include "modules/integrals.h"
#include "modules/random_draws.h"
#include "modules/driver.h"

// Get the correlation function into the integrator
CorrelationFunction * RandomDraws::corr;

STimer TotalTime;

// ====================  Computing the integrals ==================

#include "modules/compute_integral.h"
#include "modules/rescale_correlation.h"

// ================================ main() =============================


int main(int argc, char *argv[]) {

	Parameters par=Parameters(argc,argv);

    int max_no_functions=1; // required number of xi / random_draws / jackknife_weight functions
    int no_fields=1; // number of different fields used
    if(par.multi_tracers==true){
        max_no_functions=3;
        no_fields=2;
    }

    // Now read in particles to grid:
    Grid all_grid[no_fields]; // create empty grids

    for(int index=0;index<no_fields;index++){
        Float3 shift;
        Particle *orig_p;
        if (!par.make_random){
            char *filename;
            if(index==0) filename=par.fname;
            else filename=par.fname2;

            orig_p = read_particles(par.rescale, &par.np, filename, par.rstart, par.nmax);
            assert(par.np>0);

            par.perbox = compute_bounding_box(orig_p, par.np, par.rect_boxsize, par.cellsize, par.rmax, shift, par.nside);
        } else {
        // If you want to just make random particles instead:
        assert(par.np>0);
        orig_p = make_particles(par.rect_boxsize, par.np);
        par.cellsize = par.rect_boxsize.x/float(par.nside);
        // set as periodic if we make the random particles
        par.perbox = true;
        }
        if (par.qinvert) invert_weights(orig_p, par.np);
        if (par.qbalance) balance_weights(orig_p, par.np);

        // Now ready to compute!
        // Sort the particles into the grid.
        Float nofznorm=par.nofznorm;
        if(index==1) nofznorm=par.nofznorm2;
        Grid tmp_grid(orig_p, par.np, par.rect_boxsize, par.cellsize, par.nside, shift, nofznorm);

        Float grid_density = (double)par.np/tmp_grid.nf;
        printf("\n RANDOM CATALOG %d DIAGNOSTICS:\n",index+1);
        printf("Average number of particles per grid cell = %6.2f\n", grid_density);
        Float max_density = 16.0;
        if (grid_density>max_density){
            fprintf(stderr,"Average particle density exceeds maximum advised particle density (%.0f particles per cell) - exiting.\n",max_density);
            exit(1);
        }
        printf("Average number of particles per max_radius ball = %6.2f\n",
                par.np*4.0*M_PI/3.0*pow(par.rmax,3.0)/(par.rect_boxsize.x*par.rect_boxsize.y*par.rect_boxsize.z));
        if (grid_density<2){
            printf("#\n# WARNING: grid appears inefficiently fine; exiting.\n#\n");
            exit(1);
        }

        printf("# Done gridding the particles\n");
        printf("# %d particles in use, %d with positive weight\n", tmp_grid.np, tmp_grid.np_pos);
        printf("# Weights: Positive particles sum to %f\n", tmp_grid.sumw_pos);
        printf("#          Negative particles sum to %f\n", tmp_grid.sumw_neg);

        // Now save grid to global memory:
        all_grid[index].copy(&tmp_grid);

        free(orig_p); // Particles are now only stored in grid

        fflush(NULL);
    }

    // Now define all possible correlation functions and random draws:
    CorrelationFunction all_cf[max_no_functions];
    RandomDraws all_rd[max_no_functions];

    CorrelationFunction tmp_cf(par.corname,par.mbin_cf,par.mumax-par.mumin);
    all_cf[0].copy_function(&tmp_cf);
    RandomDraws tmp_rd(&tmp_cf,&par,NULL,0);
    all_rd[0].copy(&tmp_rd);

    if(par.multi_tracers==true){
        CorrelationFunction tmp_cf12(par.corname12,par.mbin_cf,par.mumax-par.mumin), tmp_cf2(par.corname2,par.mbin_cf,par.mumax-par.mumin);
        all_cf[1].copy_function(&tmp_cf2);
        all_cf[2].copy_function(&tmp_cf12);
        RandomDraws rd2(&tmp_cf2,&par,NULL,0), rd12(&tmp_cf12,&par,NULL,0);
        all_rd[1].copy(&rd2);
        all_rd[2].copy(&rd12);
    }

    // Rescale correlation functions
    rescale_correlation rescale(&par);
    rescale.refine_wrapper(&par, all_grid, all_cf, all_rd, max_no_functions);

    // Compute integrals
    compute_integral(all_grid,&par,all_cf,all_rd,1,1,1,1,1); // final digit is iteration number

    if(par.multi_tracers==true){
        compute_integral(all_grid,&par,all_cf,all_rd,1,2,1,1,2);
        compute_integral(all_grid,&par,all_cf,all_rd,1,2,2,1,3);
        compute_integral(all_grid,&par,all_cf,all_rd,1,2,1,2,4);
        compute_integral(all_grid,&par,all_cf,all_rd,1,1,2,2,5);
        compute_integral(all_grid,&par,all_cf,all_rd,2,1,2,2,6);
        compute_integral(all_grid,&par,all_cf,all_rd,2,2,2,2,7);
    }

    rusage ru;
    getrusage(RUSAGE_SELF, &ru);

    fprintf(stderr,"# user CPU time used: %ld \n# system CPU time used: %ld \n# maximum resident set size: %ld \n# integral shared memory size: %ld \n# integral unshared data size: %ld \n# integral unshared stack size: %ld \n# page reclaims (soft page faults): %ld \n# page faults (hard page faults): %ld \n# swaps: %ld \n# block input operations: %ld \n# block output operations: %ld \n# IPC messages sent: %ld \n# IPC messages received: %ld \n# signals received: %ld \n# voluntary context switches: %ld \n# involuntary context switches: %ld \n",ru.ru_utime.tv_sec,ru.ru_stime.tv_sec,ru.ru_maxrss,ru.ru_ixrss,ru.ru_idrss,ru.ru_isrss,ru.ru_minflt,ru.ru_majflt,ru.ru_nswap,ru.ru_inblock,ru.ru_oublock,ru.ru_msgsnd,ru.ru_msgrcv,ru.ru_nsignals,ru.ru_nvcsw,ru.ru_nivcsw);

    return 0;
}
