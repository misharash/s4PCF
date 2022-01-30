#ifndef PAIRS_H
#define PAIRS_H

// ========================== Accumulate the two-pcf pair counts
// ================

class Pairs {
   public:
    double* xi0;

   private:
    double* xi2;

   private:
    double empty[8];  // Just to try to keep the threads from working on similar
                      // memory

   public:
    Pairs() {
        // Initialize the binning
        int ec = 0;
        ec += posix_memalign((void**)&xi0, PAGE, sizeof(double) * NBIN_SHORT);
        ec += posix_memalign((void**)&xi2, PAGE, sizeof(double) * NBIN_SHORT);
        assert(ec == 0);
        for (int j = 0; j < NBIN_SHORT; j++) {
            xi0[j] = 0;
            xi2[j] = 0;
        }
        empty[0] = 0.0;  // To avoid a warning
    }
    ~Pairs() {
        free(xi0);
        free(xi2);
    }

    inline void load(Float* xi0ptr, Float* xi2ptr) {
        for (int j = 0; j < NBIN_SHORT; j++) {
            xi0[j] = xi0ptr[j];
            xi2[j] = xi2ptr[j];
        }
    }

    inline void save(Float* xi0ptr, Float* xi2ptr) {
        for (int j = 0; j < NBIN_SHORT; j++) {
            xi0ptr[j] = xi0[j];
            xi2ptr[j] = xi2[j];
        }
    }

    inline void add(int b, Float dz, Float w) {
        // Add up the weighted pair for the monopole and quadrupole correlation
        // function
        xi0[b] += w;
        xi2[b] += w * (3.0 * dz * dz - 1) * 0.5;
    }

    void sum_power(Pairs* p) {
        // Just add up all of the threaded pairs into the zeroth element
        for (int i = 0; i < NBIN_SHORT; i++) {
            xi0[i] += p->xi0[i];
            xi2[i] += p->xi2[i];
        }
    }
    //
    //   void report_pairs() {
    //   for (int j=0; j<NBIN_SHORT; j++) {
    //     printf("Pairs %2d %9.0f %9.0f\n",
    // 		j, xi0[j], xi2[j]);
    // }
    //   }

    void save_pairs(char* out_string, Float rmin_short, Float rmax_short, Float sumw) {
        // Print the output isotropic 2PCF counts to file

        // Create output directory if not in existence
        const char* out_dir;
        out_dir = "output";
        if (mkdir(out_dir, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH) == 0) {
            printf("\nCreating output directory\n");
        }

        // First create output files
        char out_name[1000];
        snprintf(out_name, sizeof out_name, "output/%s_2pcf.txt", out_string);
        FILE* OutFile = fopen(out_name, "w");

        // Print some useful information
        fprintf(OutFile, "## Bins: %d\n", NBIN_SHORT);
        fprintf(OutFile, "## Minimum Radius = %.2e\n", rmin_short);
        fprintf(OutFile, "## Maximum Radius = %.2e\n", rmax_short);
        fprintf(OutFile, "## Format: Row 1 = radial bin 1, Row 2 = xi^a\n");

        // First print the indices of the first radial bin
        for (int i = 0; i < NBIN_SHORT; i++)
            fprintf(OutFile, "%2d\t", i);
        fprintf(OutFile, " \n");

        // Now print the 2PCF
        Float norm = pow(sumw, -2); // normalize by sum of (positive) weights squared
        for (int i = 0; i < NBIN_SHORT; i++)
            fprintf(OutFile, "%le\t", xi0[i]*norm);
        fprintf(OutFile, "\n");

        fflush(NULL);

        // Close open files
        fclose(OutFile);

        printf("\n2PCF Output saved to %s\n", out_name);
    }
};

#endif