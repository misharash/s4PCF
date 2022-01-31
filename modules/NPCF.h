#ifndef NPCF_H
#define NPCF_H

// ========================== Here are all of the cross-powers ============

class NPCF {
    // This should accumulate the NPCF contributions, for all combination of
    // bins.
   public:

// Array to hold the 3PCF
#define N3PCF (NBIN_SHORT * (NBIN_SHORT + 1) / 2)  // only compute bin1 <= bin2
    Float threepcf[N3PCF];

// Sizes of 4pcf array
#define N4PCF NBIN_LONG* NBIN_SHORT*(NBIN_SHORT + 1) / 2

    // Array to hold the 4PCF
    Float fourpcf[N4PCF];

    int N4PCF_used = 0; // number of 4PCF actually in use, can't be larger than N4PCF, will be incremented later

    int fourpcf_bins[3][N4PCF]; // 4PCF index to 3 bin numbers mapping

    int fourpcf_bin_number[NBIN_LONG][NBIN_SHORT][NBIN_SHORT]; // 3 bin numbers to 4PCF index mapping, uses -1 when bin is illegal

    void reset() {
        // Zero out the array on construction.

        for (int x = 0; x < N3PCF; x++)
            threepcf[x] = 0.0;

        // Initialize array to zero
        for (int x = 0; x < N4PCF; x++)
            fourpcf[x] = 0.0;

        return;
    }

    void calc_4pcf_indices(Float rmin_short, Float rmax_short, Float rmin_long, Float rmax_long) {
        // Set up index mapping
        // Check the triangle condition if needed
        double eps = 1e-6; // some small tolerance for checking the condition
        for (int i = 0; i < NBIN_LONG; i++) {
            Float ri_min = rmin_long + i*(rmax_long-rmin_long)/NBIN_LONG; // minimal long radius in i'th bin
            for (int j = 0; j < NBIN_SHORT; j++) {
                Float rj_max = rmin_short + (j+1)*(rmax_short-rmin_short)/NBIN_SHORT; // maximal radius in j'th bin
                for (int k = j; k < NBIN_SHORT; k++) { // enough to do k>=j
                    Float rk_max = rmin_short + (k+1)*(rmax_short-rmin_short)/NBIN_SHORT; // maximal radius in k'th bin
                    if (PREVENT_TRIANGLES && (rj_max + rk_max > ri_min + eps)) { // triangle is possible and we prevent it
                        fourpcf_bin_number[i][j][k] = fourpcf_bin_number[i][k][j] = -1; // bin is illegal, set that symmetrically
                    }
                    else {
                        fourpcf_bin_number[i][j][k] = fourpcf_bin_number[i][k][j] = N4PCF_used; // set bin number symmetrically by last two numbers
                        fourpcf_bins[0][N4PCF_used] = i;
                        fourpcf_bins[1][N4PCF_used] = j;
                        fourpcf_bins[2][N4PCF_used] = k;
                        N4PCF_used++; // increment number of used bins
                    }
                }
            }
        }
        assert(N4PCF_used <= N4PCF);
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

    void save_power(char* out_string,
                    Float rmin_short,
                    Float rmax_short,
                    Float rmin_long,
                    Float rmax_long,
                    Float sumw) {
        // Print the output NPCF counts to file

        // SAVE 3PCF

        // First create output files
        char out_name[1000];
        snprintf(out_name, sizeof out_name, "output/%s_3pcf.txt", out_string);
        FILE* OutFile = fopen(out_name, "w");

        // Print some useful information
        fprintf(OutFile, "## Bins: %d\n", NBIN_SHORT);
        fprintf(OutFile, "## Minimum Radius = %.2e\n", rmin_short);
        fprintf(OutFile, "## Maximum Radius = %.2e\n", rmax_short);
        fprintf(OutFile,
                "## Format: Row 1 = radial bin 1, Row 2 = radial bin 2, Row "
                "3 = zeta_ell^ab\n");

        // First print the indices of the first radial bin
        for (int i = 0; i < NBIN_SHORT; i++) {
            for (int j = i; j < NBIN_SHORT; j++)
                fprintf(OutFile, "%2d\t", i);
        }
        fprintf(OutFile, " \n");

        // Print the indices of the second radial bin
        for (int i = 0; i < NBIN_SHORT; i++) {
            for (int j = i; j < NBIN_SHORT; j++)
                fprintf(OutFile, "%2d\t", j);
        }
        fprintf(OutFile, "\n");

        // Now print the 3PCF
        Float norm = pow(sumw, -3); // normalize by sum of (positive) weights cubed
        for (int ct = 0; ct < N3PCF; ct++)
            fprintf(OutFile, "%le\t", threepcf[ct]*norm);
        
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
            fprintf(OutFile2, "## Short side bins: %d\n", NBIN_SHORT);
            fprintf(OutFile2, "## Long side bins: %d\n", NBIN_LONG);
            fprintf(OutFile2, "## Minimum Radius for short side = %.2e\n", rmin_short);
            fprintf(OutFile2, "## Maximum Radius for short side = %.2e\n", rmax_short);
            fprintf(OutFile2, "## Minimum Radius for long side = %.2e\n",
                    rmin_long);
            fprintf(OutFile2, "## Maximum Radius for long side = %.2e\n",
                    rmax_long);
            fprintf(OutFile2,
                    "## Format: Row 1 = radial bin 1 (long side), Row 2 = radial bin 2 (short side), "
                    "Row 3 = radial bin 3 (short side), Row 4 = zeta_l1l2l3^abc\n");

            // First print the indices of the radial bins
            for (int i = 0; i < N4PCF_used; i++)
                fprintf(OutFile2, "%2d\t", fourpcf_bins[0][i]);
            fprintf(OutFile2, "\n");

            for (int i = 0; i < N4PCF_used; i++)
                fprintf(OutFile2, "%2d\t", fourpcf_bins[1][i]);
            fprintf(OutFile2, "\n");

            for (int i = 0; i < N4PCF_used; i++)
                fprintf(OutFile2, "%2d\t", fourpcf_bins[2][i]);
            fprintf(OutFile2, "\n");

            // Now print the 4PCF.
            Float norm = pow(sumw, -4); // normalize by sum of (positive) weights to the 4th power
            for (int i = 0; i < N4PCF_used; i++)
                fprintf(OutFile2, "%le\t", fourpcf[i]*norm);
            fprintf(OutFile2, "\n");
            fflush(NULL);

            // Close open files
            fclose(OutFile2);

            printf("4PCF Output saved to %s\n", out_name2);
        }
    }

    inline int getbin_pair(int i, int j) { // needs i <= j
        return j + NBIN_SHORT*i - i * (i-1) / 2;
    }

    inline void excl_3pcf(int bin, Float wprod) {
        // wprod is product of weights

        // delete 3PCF self-count

        int index = getbin_pair(bin, bin);
        threepcf[index] -= wprod;

        return;
    }

    inline void add_3pcf(Pairs* pairs, Float wp) {
        // wp is the primary galaxy weight

        // COMPUTE 3PCF CONTRIBUTIONS

        for (int i = 0, ct = 0; i < NBIN_SHORT; i++) {
            for (int j = i; j < NBIN_SHORT; j++, ct++) {
                threepcf[ct] = pairs->xi0[i] * pairs->xi0[j] / wp;
            }
        }

        return;
    }

#if (!PREVENT_TRIANGLES)
    inline void excl_4pcf_triangle(int bin_long, int bin, int bin2, Float wprod) {
        // COMPUTE 4PCF CONTRIBUTIONS

        if (bin <= bin2)
            fourpcf[fourpcf_bin_number[bin_long][bin][bin2]] -= wprod;

        return;
    }

    inline void excl_4pcf_doubleside(Pairs* pair, int bin_long, int bin, Float wprod) {
        // COMPUTE 4PCF CONTRIBUTIONS

        // Iterate over second bin
        for (int bin2 = 0; bin2 <= bin; bin2++)
            fourpcf[fourpcf_bin_number[bin_long][bin][bin2]] -= pair->xi0[bin2] * wprod;
        for (int bin2 = bin; bin2 < NBIN_SHORT; bin2++) // go over bin=bin2 second time for consistency
            fourpcf[fourpcf_bin_number[bin_long][bin][bin2]] -= pair->xi0[bin2] * wprod;
        // End of radial binning loops

        return;
    }

    inline void excl_4pcf_tripleside(int bin_long, int bin, Float wprod) {
        // COMPUTE 4PCF CONTRIBUTIONS

        fourpcf[fourpcf_bin_number[bin_long][bin][bin]] += wprod * wprod;

        return;
    }
#endif

    inline void add_4pcf(Pairs* pair1, Pairs* pair2, int bin_long) {
        // COMPUTE 4PCF CONTRIBUTIONS

        // Iterate over second bin
        for (int j = 0; j < NBIN_SHORT; j++) {
            // Iterate over final bin and advance the 4PCF array counter
            for (int k = j; k < NBIN_SHORT; k++) {
                int bin_number = fourpcf_bin_number[bin_long][j][k];
                #if (PREVENT_TRIANGLES)
                if (bin_number < 0) break; // skip illegal bin, then all next ones are illegal too. No illegal bins if no triangle prevention.
                #endif
                fourpcf[bin_number] += pair1->xi0[j] * pair2->xi0[k] + pair1->xi0[k] * pair2->xi0[j];
            }
        }
        // End of radial binning loops

        return;
    }

};  // end NPCF class

#endif
