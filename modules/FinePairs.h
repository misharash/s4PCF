#ifndef FINE_PAIRS_H
#define FINE_PAIRS_H

// Accumulate the anisotropic, finer two-pcf pair counts

class FinePairs {
   public:
    double* xi;

   private:
    double empty[8];  // Just to try to keep the threads from working on similar
                      // memory

   public:
    FinePairs() {
        // Initialize the binning
        int ec = 0;
        ec += posix_memalign((void**)&xi, PAGE, sizeof(double) * NBIN_CF * MBIN_CF);
        assert(ec == 0);
        for (int j = 0; j < NBIN_CF*MBIN_CF; j++)
            xi[j] = 0;
        empty[0] = 0.0;  // To avoid a warning
    }
    ~FinePairs() {
        free(xi);
    }

    inline void load(Float* xiptr) {
        for (int j = 0; j < NBIN_CF*MBIN_CF; j++)
            xi[j] = xiptr[j];
    }

    inline void save(Float* xiptr) {
        for (int j = 0; j < NBIN_CF*MBIN_CF; j++)
            xiptr[j] = xi[j];
    }

    inline void add(int b, Float dz, Float w) {
        // Add up the weighted pair for the fine anisotropic correlation function
        int mubin = floor(fabs(dz)*MBIN_CF);
        xi[b * MBIN_CF + mubin] += w;
    }

    void sum_power(FinePairs* p) {
        // Just add up all of the threaded pairs into the zeroth element
        for (int j = 0; j < NBIN_CF*MBIN_CF; j++)
            xi[j] += p->xi[j];
    }

    void save_pairs(char* out_string, Float rmin_cf, Float rmax_cf, Float sumw) {
        // Print the output anisotropic 2PCF counts to file

        // First create output files
        char out_name[1000];
        snprintf(out_name, sizeof out_name, "output/%s_fine2pcf.txt", out_string);
        FILE* OutFile = fopen(out_name, "w");

        // First print the midpoints of radial bins
        for (int i = 0; i < NBIN_CF; i++)
            fprintf(OutFile, "%2e\t", rmin_cf + (rmax_cf - rmin_cf) * ((Float)i + 0.5) / NBIN_CF);
        fprintf(OutFile, "\n");

        // First print the midpoints of mu bins
        for (int i = 0; i < MBIN_CF; i++)
            fprintf(OutFile, "%2e\t", ((Float)i + 0.5) / MBIN_CF);
        fprintf(OutFile, "\n");

        // Now print the 2PCF
        Float norm = pow(sumw, -2); // normalize by sum of (positive) weights squared
        for (int i = 0; i < NBIN_CF; i++) {
            for (int j = 0; j < MBIN_CF; j++)
                fprintf(OutFile, "%le\t", xi[i * MBIN_CF + j]*norm);
            fprintf(OutFile, "\n");
        }

        fflush(NULL);

        // Close open files
        fclose(OutFile);

        printf("\n2PCF Output saved to %s\n", out_name);
    }
};

#endif