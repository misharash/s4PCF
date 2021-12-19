#ifndef NPCF_H
#define NPCF_H

// ========================== Here are all of the cross-powers ============

class NPCF {
    // This should accumulate the NPCF contributions, for all combination of
    // bins.
   public:
    STimer AddTimer3, ExclTimer3;

// Array to hold the 3PCF
#define N3PCF (NBIN * (NBIN + 1) / 2)  // only compute bin1 <= bin2
    Float threepcf[N3PCF];

    STimer AddTimer4, ExclTimer4_doubleside, ExclTimer4_tripleside, ExclTimer4_triangle;

// Sizes of 4pcf array
#define N4PCF NBIN_LONG* NBIN*(NBIN + 1) / 2

    // Array to hold the 4PCF
    Float fourpcf[N4PCF];

    void reset() {
        // Zero out the array on construction.

        for (int x = 0; x < N3PCF; x++)
            threepcf[x] = 0.0;

        // Initialize array to zero
        for (int x = 0; x < N4PCF; x++)
            fourpcf[x] = 0.0;

        return;
    }

    NPCF() {
        reset();
        return;
    }
    ~NPCF() {}

    void sum_power(NPCF* c) {
        // Just add up all of the threaded power into the zeroth element

        for (int x = 0; x < N3PCF; x++)
            threepcf[x] += c->threepcf[x];

        for (int x = 0; x < N4PCF; x++)
            fourpcf[x] += c->fourpcf[x];
    }

    void report_timings() {
        /// Report the NPCF timing measurements (for a single CPU).

        printf("\n# Single CPU Timings");
        printf("\n3PCF self-counts exclusion: %.3f s", ExclTimer3.Elapsed());
        printf("\n3PCF addition: %.3f s", AddTimer3.Elapsed());
        printf("\n4PCF triangle self-counts exclusion: %.3f s", ExclTimer4_triangle.Elapsed());
        printf("\n4PCF double-side self-counts exclusion: %.3f s", ExclTimer4_doubleside.Elapsed());
        printf("\n4PCF triple-side self-counts exclusion: %.3f s", ExclTimer4_tripleside.Elapsed());
        printf("\n4PCF addition: %.3f s", AddTimer4.Elapsed());
        printf("\n");
    }

    void save_power(char* out_string,
                    Float rmin,
                    Float rmax,
                    Float rmin_long,
                    Float rmax_long) {
        // Print the output NPCF counts to file

        // SAVE 3PCF

        // First create output files
        char out_name[1000];
        snprintf(out_name, sizeof out_name, "output/%s_3pcf.txt", out_string);
        FILE* OutFile = fopen(out_name, "w");

        // Print some useful information
        fprintf(OutFile, "## Bins: %d\n", NBIN);
        fprintf(OutFile, "## Minimum Radius = %.2e\n", rmin);
        fprintf(OutFile, "## Maximum Radius = %.2e\n", rmax);
        fprintf(OutFile,
                "## Format: Row 1 = radial bin 1, Row 2 = radial bin 2, Row "
                "3 = zeta_ell^ab\n");

        // First print the indices of the first radial bin
        for (int i = 0; i < NBIN; i++) {
            for (int j = i; j < NBIN; j++)
                fprintf(OutFile, "%2d\t", i);
        }
        fprintf(OutFile, " \n");

        // Print the indices of the second radial bin
        for (int i = 0; i < NBIN; i++) {
            for (int j = i; j < NBIN; j++)
                fprintf(OutFile, "%2d\t", j);
        }
        fprintf(OutFile, "\n");

        // Now print the 3PCF
        for (int ct = 0; ct < N3PCF; ct++)
            fprintf(OutFile, "%le\t", threepcf[ct]);
        
        fprintf(OutFile, "\n");

        fflush(NULL);

        // Close open files
        fclose(OutFile);

        printf("3PCF Output saved to %s\n", out_name);

        {
            // SAVE 4PCF

            // First create output files
            char out_name2[1000];
            snprintf(out_name2, sizeof out_name2, "output/%s_4pcf.txt",
                     out_string);
            FILE* OutFile2 = fopen(out_name2, "w");

            // Print some useful information
            fprintf(OutFile2, "## Bins: %d\n", NBIN);
            fprintf(OutFile2, "## Minimum Radius = %.2e\n", rmin);
            fprintf(OutFile2, "## Maximum Radius = %.2e\n", rmax);
            fprintf(OutFile2, "## Minimum Radius for long side = %.2e\n",
                    rmin_long);
            fprintf(OutFile2, "## Maximum Radius for long side = %.2e\n",
                    rmax_long);
            fprintf(OutFile2,
                    "## Format: Row 1 = radial bin 1 (long side), Row 2 = radial bin 2, "
                    "Row 3 = radial bin 3, Row 4 = zeta_l1l2l3^abc\n");

            // First print the indices of the radial bins
            for (int i = 0; i < NBIN_LONG; i++) {
                for (int j = 0; j < NBIN; j++) {
                    for (int k = j; k < NBIN; k++) {
                        fprintf(OutFile2, "%2d\t", i);
                    }
                }
            }
            fprintf(OutFile2, "\n");

            for (int i = 0; i < NBIN_LONG; i++) {
                for (int j = 0; j < NBIN; j++) {
                    for (int k = j; k < NBIN; k++) {
                        fprintf(OutFile2, "%2d\t", j);
                    }
                }
            }
            fprintf(OutFile2, "\n");

            for (int i = 0; i < NBIN_LONG; i++) {
                for (int j = 0; j < NBIN; j++) {
                    for (int k = j; k < NBIN; k++) {
                        fprintf(OutFile2, "%2d\t", k);
                    }
                }
            }
            fprintf(OutFile2, "\n");

            // Now print the 4PCF.
            for (int i = 0; i < N4PCF; i++)
                fprintf(OutFile2, "%le\t", fourpcf[i]);
            fprintf(OutFile2, "\n");
            fflush(NULL);

            // Close open files
            fclose(OutFile2);

            printf("4PCF Output saved to %s\n", out_name2);
        }
    }

    inline void excl_3pcf(int bin, Float wprod) {
        // wprod is product of weights

        // delete 3PCF self-count
        ExclTimer3.Start();

        int i = bin, j = bin;
        int index = j + NBIN*i - i * (i-1) / 2;
        threepcf[index] -= wprod;

        ExclTimer3.Stop();

        return;
    }

    inline void add_3pcf(Pairs* pairs, Float wp) {
        // wp is the primary galaxy weight

        // COMPUTE 3PCF CONTRIBUTIONS
        AddTimer3.Start();

        for (int i = 0, ct = 0; i < NBIN; i++) {
            for (int j = i; j < NBIN; j++, ct++) {
                threepcf[ct] = pairs->xi0[i] * pairs->xi0[j] / wp;
            }
        }

        AddTimer3.Stop();

        return;
    }

    inline void excl_4pcf_triangle(int bin_long, int bin, int bin2, Float wprod) {
        // COMPUTE 4PCF CONTRIBUTIONS

        ExclTimer4_triangle.Start();

        int i = bin, j = bin2; // i <= j
        if (i <= j) {
            int index = j + NBIN*i - i * (i-1) / 2;
            fourpcf[bin_long * N3PCF + index] -= wprod;
        }

        ExclTimer4_triangle.Stop();

        return;
    }

    inline void excl_4pcf_doubleside(Pairs* pair, int bin_long, int bin, Float wprod) {
        // COMPUTE 4PCF CONTRIBUTIONS

        ExclTimer4_doubleside.Start();

        // Iterate over second bin
        for (int bin2 = 0; bin2 <= bin; bin2++) {
            int i = bin2, j = bin; // i <= j
            int index = j + NBIN*i - i * (i-1) / 2;
            fourpcf[bin_long * N3PCF + index] -= pair->xi0[bin2] * wprod;
        }
        for (int bin2 = bin; bin2 < NBIN; bin2++) { // go over bin=bin2 second time for consistency
            int i = bin, j = bin2; // i <= j
            int index = j + NBIN*i - i * (i-1) / 2;
            fourpcf[bin_long * N3PCF + index] -= pair->xi0[bin2] * wprod;
        }
        // End of radial binning loops
        ExclTimer4_doubleside.Stop();

        return;
    }

    inline void excl_4pcf_tripleside(int bin_long, int bin, Float wprod) {
        // COMPUTE 4PCF CONTRIBUTIONS

        ExclTimer4_tripleside.Start();

        int i = bin, j = bin; // i <= j
        int index = j + NBIN*i - i * (i-1) / 2;
        fourpcf[bin_long * N3PCF + index] += wprod * wprod;
        
        ExclTimer4_tripleside.Stop();

        return;
    }

    inline void add_4pcf(Pairs* pair1, Pairs* pair2, int bin_long) {
        // COMPUTE 4PCF CONTRIBUTIONS

        AddTimer4.Start();

        // Iterate over second bin
        for (int j = 0, bin_index = 0; j < NBIN; j++) {
            // Iterate over final bin and advance the 4PCF array counter
            for (int k = j; k < NBIN; k++) {
                fourpcf[bin_long * N3PCF + bin_index++] +=
                    pair1->xi0[j] * pair2->xi0[k] +
                    pair1->xi0[k] * pair2->xi0[j];
            }
        }
        // End of radial binning loops
        AddTimer4.Stop();

        return;
    }

};  // end NPCF class

#endif
