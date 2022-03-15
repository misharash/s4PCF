#ifndef COMPUTE_PAIRS_H
#define COMPUTE_PAIRS_H

// ====================  Computing the pairs ==================

void compute_pairs(Grid* grid,
                   Float rmin_short,
                   Float rmax_short,
                   Float rmin_long,
                   Float rmax_long,
                   Float rmin_cf,
                   Float rmax_cf,
                   int np) {
    int maxsep_short = ceil(rmax_short / grid->cellsize);  // Maximum distance we must search for short side
    int maxsep_long_or_cf = ceil(fmax(rmax_long, rmax_cf) / grid->cellsize);  // Maximum distance we must search for long side
    int ne;
    Float rmax_short2 = rmax_short * rmax_short;
    Float rmin_short2 = rmin_short * rmin_short;  // rmax2*1e-12;    // Just an underflow guard
    Float rmax_long2 = rmax_long * rmax_long;
    Float rmin_long2 = rmin_long * rmin_long;
    Float rmax_cf2 = rmax_cf * rmax_cf;
    Float rmin_cf2 = rmin_cf * rmin_cf;
    uint64 cnt = 0, cnt2 = 0, cnt3 = 0;

    Pairs* pairs_i = new Pairs[np];

    // Easy to multi-thread this top loop!
    // But some cells have trivial amounts of work, so we will first make a list
    // of the work. Including the empty cells appears to fool the dynamic thread
    // allocation sometimes.

    STimer accpairs,
        powertime;  // measure the time spent accumulating powers for multipoles
    // We're going to loop only over the non-empty cells.

    accpairs.Start();

#ifdef OPENMP
#pragma omp parallel for schedule(dynamic, 8) reduction(+ : cnt)
#endif

    for (ne = 0; ne < grid->nf; ne++) {
        int n = grid->filled[ne];  // Fetch the cell number

        // Decide which thread we are in
#ifdef OPENMP
        int thread = omp_get_thread_num();
        assert(omp_get_num_threads() <= MAXTHREAD);
        if (ne == 0)
            printf("# Running on %d threads.\n", omp_get_num_threads());
#else
        int thread = 0;
        if (ne == 0)
            printf("# Running single threaded.\n");
#endif
        if (int(ne % 1000) == 0)
            printf("Computing cell %d of %d on thread %d (1st loop)\n", ne, grid->nf, thread);

        // Loop over primary cells.
        Cell primary = grid->c[n];
        integer3 prim_id = grid->cell_id_from_1d(n);

        // continue; // To skip all of the list-building and summations.
        // Everything else takes negligible time
        // Now we need to loop over all primary particles in this cell
        for (int j = primary.start; j < primary.start + primary.np; j++) {
            Float primary_w = grid->p[j].w;

            #if (PERIODIC)
            // do analytic pair counts for random secondary particles, i.e. take limit of their number to infinity
            for (int bin = 0; bin < NBIN_SHORT; bin++) {
                Float rmin_bin = rmin_short + bin * (rmax_short - rmin_short) / NBIN_SHORT;
                Float rmax_bin = rmin_short + (bin + 1) * (rmax_short - rmin_short) / NBIN_SHORT;
                pairs_i[j].add(bin, primary_w * (grid->sumw_neg) / pow(grid->max_boxsize, 3) * 4 * M_PI / 3 * (pow(rmax_bin, 3) - pow(rmin_bin, 3)));
            }
            #endif

            // Then loop over secondaries, cell-by-cell
            integer3 delta;
            for (delta.x = -maxsep_short; delta.x <= maxsep_short; delta.x++)
                for (delta.y = -maxsep_short; delta.y <= maxsep_short; delta.y++)
                    for (delta.z = -maxsep_short; delta.z <= maxsep_short; delta.z++) {
                        const int samecell = (delta.x == 0 && delta.y == 0 && delta.z == 0) ? 1 : 0;

                        // Check that the cell is in the grid!
                        int tmp_test = grid->test_cell(prim_id + delta);
                        if (tmp_test < 0)
                            continue;
                        Cell sec = grid->c[tmp_test];

                        // Define primary position
                        Float3 ppos = grid->p[j].pos;
                        #if PERIODIC
                        ppos -= grid->cell_sep(delta);
                        #endif

                        // This is the position of the particle as viewed
                        // from the secondary cell. Now loop over the
                        // particles in this secondary cell
                        for (int k = sec.start; k < sec.start + sec.np;
                                k++) {
                            // Now we're considering these two particles!
                            if (samecell && j == k)
                                continue;  // Exclude self-count
                            Float3 dx = grid->p[k].pos - ppos;
                            Float norm2 = dx.norm2();
                            // Check if this is in the correct binning
                            // ranges
                            if (norm2 < rmax_short2 && norm2 > rmin_short2)
                                cnt++;
                            else
                                continue;

                            // Now what do we want to do with the pair?
                            norm2 = sqrt(norm2);  // Now just radius
                            // Find the radial bin
                            int bin = floor((norm2 - rmin_short) / (rmax_short - rmin_short) * NBIN_SHORT);

                            // Accumulate the 2-pt correlation function
                            // We include the weight for each pair
                            #if (PERIODIC)
                            if (grid->p[k].w > 0) // in periodic box, count only secondary data particles, randoms added earlier analytically
                            #endif
                            pairs_i[j].add(bin, grid->p[k].w * primary_w);
                            // Exclude self-counts from 3PCF
                            npcf[thread].excl_3pcf(bin, grid->p[k].w * grid->p[k].w * primary_w);

                            // Exclude triangular self-counts from 4PCF
                            #if (!PREVENT_TRIANGLES && !IGNORE_TRIANGLES)
                            if (rmin_long < 2*rmax_short) {
                                integer3 delta2;
                                for (delta2.x = -maxsep_short; delta2.x <= maxsep_short; delta2.x++)
                                    for (delta2.y = -maxsep_short; delta2.y <= maxsep_short; delta2.y++)
                                        for (delta.z = -maxsep_short; delta.z <= maxsep_short; delta2.z++) {
                                            // Check that the cell is in the grid!
                                            int tmp_test = grid->test_cell(prim_id + delta2);
                                            if (tmp_test < 0)
                                                continue;
                                            Cell third = grid->c[tmp_test];

                                            // Define primary position
                                            Float3 ppos2 = grid->p[j].pos;
                                            #if PERIODIC
                                            ppos2 -= grid->cell_sep(delta2);
                                            #endif

                                            // Define secondary position
                                            Float3 spos2 = grid->p[k].pos;
                                            #if PERIODIC
                                            spos2 -= grid->cell_sep(delta2) - grid->cell_sep(delta);
                                            #endif
                                            // This is the position of the particle as viewed
                                            // from the secondary cell. Now loop over the
                                            // particles in this secondary cell
                                            for (int l = third.start; l < third.start + third.np;
                                                l++) {
                                                // Now we're considering these two particles!
                                                if ((j==l) || (k==l))
                                                    continue;  // Exclude self-count
                                                Float3 dx = grid->p[l].pos - ppos2;
                                                Float norm_2 = dx.norm2();
                                                Float3 dx_l = grid->p[l].pos - spos2;
                                                Float norm_l2 = dx_l.norm2();
                                                // Check if this is in the correct binning
                                                // ranges
                                                if (norm_2 < rmax_short2 && norm_2 > rmin_short2 && norm_l2 < rmax_long2 && norm_l2 > rmin_long2)
                                                    cnt3++;
                                                else
                                                    continue;

                                                // Now what do we want to do with the pair?
                                                norm_2 = sqrt(norm_2);  // Now just radius
                                                norm_l2 = sqrt(norm_l2);
                                                // Find the radial bins
                                                int bin2 = floor((norm_2 - rmin_short) / (rmax_short - rmin_short) * NBIN_SHORT);
                                                // if (bin2 < bin) continue; // count triples only in one order
                                                int bin_long = floor((norm_l2 - rmin_long) / (rmax_long - rmin_long) * NBIN_LONG);

                                                // Exclude self-counts from 4PCF
                                                npcf[thread].excl_4pcf_triangle(bin_long, bin, bin2,
                                                                                grid->p[l].w * grid->p[k].w * primary_w * primary_w);
                                            }  // Done with this secondary particle
                                        }      // Done with this delta.z loop
                                // done with delta.y loop
                                // done with delta.x loop
                            }
                            #endif

                            // Exclude triple-side self-counts from 4PCF
                            #if (!PREVENT_TRIANGLES && !IGNORE_TRIANGLES)
                            if (rmin_long < rmax_short) {
                                if (norm2 >= rmax_long2 || norm2 <= rmin_long2)
                                    continue;
                                int bin_long = floor((norm2 - rmin_long) / (rmax_long - rmin_long) * NBIN_LONG);
                                npcf[thread].excl_4pcf_tripleside(bin_long, bin, grid->p[k].w * primary_w);
                            }
                            #endif
                        }  // Done with this secondary particle
                    }      // Done with this delta.z loop
            // done with delta.y loop
            // done with delta.x loop
            pairs[thread].sum_power(pairs_i + j);

            // Now exclude 4pcf double-side self-counts
            #if (!PREVENT_TRIANGLES && !IGNORE_TRIANGLES)
            if (rmin_long < rmax_short) {
                for (delta.x = -maxsep_short; delta.x <= maxsep_short; delta.x++)
                    for (delta.y = -maxsep_short; delta.y <= maxsep_short; delta.y++)
                        for (delta.z = -maxsep_short; delta.z <= maxsep_short; delta.z++) {
                            const int samecell =
                                (delta.x == 0 && delta.y == 0 && delta.z == 0)
                                    ? 1
                                    : 0;

                            // Check that the cell is in the grid!
                            int tmp_test = grid->test_cell(prim_id + delta);
                            if (tmp_test < 0)
                                continue;
                            Cell sec = grid->c[tmp_test];

                            // Define primary position
                            Float3 ppos = grid->p[j].pos;
                            #if PERIODIC
                            ppos -= grid->cell_sep(delta);
                            #endif

                            // This is the position of the particle as viewed
                            // from the secondary cell. Now loop over the
                            // particles in this secondary cell
                            for (int k = sec.start; k < sec.start + sec.np;
                                k++) {
                                // Now we're considering these two particles!
                                if (samecell && j == k)
                                    continue;  // Exclude self-count
                                Float3 dx = grid->p[k].pos - ppos;
                                Float norm2 = dx.norm2();
                                // Check if this is in the correct binning
                                // ranges
                                if (norm2 >= rmax_short2 || norm2 <= rmin_short2)
                                    continue;

                                // Now what do we want to do with the pair?
                                norm2 = sqrt(norm2);  // Now just radius
                                // Find the radial bin
                                int bin = floor((norm2 - rmin_short) / (rmax_short - rmin_short) * NBIN_SHORT);

                                // Define x/r,y/r,z/r
                                dx = dx / norm2;

                                // continue;   // Skip pairs and multipoles

                                // Exclude self-counts from 4PCF
                                if (norm2 >= rmax_long2 || norm2 <= rmin_long2)
                                    continue;
                                // Find the radial bin for long side
                                int bin_long = floor((norm2 - rmin_long) / (rmax_long - rmin_long) * NBIN_LONG);
                                npcf[thread].excl_4pcf_doubleside(pairs_i + j, bin_long, bin, grid->p[k].w * primary_w);
                            }  // Done with this secondary particle
                        }      // Done with this delta.z loop
                // done with delta.y loop
                // done with delta.x loop
            }
            #endif

            // Now combine pair counts into 3pcf counts
            npcf[thread].add_3pcf(pairs_i + j, primary_w);
        }  // Done with this primary particle

    }  // Done with this primary cell (first loop), end of omp pragma

    accpairs.Stop();

    // Second loop, to make sure pair counts are ready

    powertime.Start();

#ifdef OPENMP
#pragma omp parallel for schedule(dynamic, 8) reduction(+ : cnt)
#endif

    for (ne = 0; ne < grid->nf; ne++) {
        int n = grid->filled[ne];  // Fetch the cell number

        // Decide which thread we are in
#ifdef OPENMP
        int thread = omp_get_thread_num();
        assert(omp_get_num_threads() <= MAXTHREAD);
#else
        int thread = 0;
#endif
        if (int(ne % 1000) == 0)
            printf("Computing cell %d of %d on thread %d (2nd loop)\n", ne, grid->nf,
                   thread);

        // Loop over primary cells.
        Cell primary = grid->c[n];
        integer3 prim_id = grid->cell_id_from_1d(n);

        // continue; // To skip all of the list-building and summations.
        // Everything else takes negligible time
        // Now we need to loop over all primary particles in this cell
        for (int j = primary.start; j < primary.start + primary.np; j++) {
            Float primary_w = grid->p[j].w;
            // Then loop over secondaries, cell-by-cell
            integer3 delta;

            // Now combine pair counts into 4pcf counts

            // This is done on CPU - calculate add_to_power here
            // Acumulate powers here - code in NPCF.h
            for (delta.x = -maxsep_long_or_cf; delta.x <= maxsep_long_or_cf;
                    delta.x++)
                for (delta.y = -maxsep_long_or_cf; delta.y <= maxsep_long_or_cf;
                        delta.y++)
                    for (delta.z = -maxsep_long_or_cf; delta.z <= maxsep_long_or_cf;
                            delta.z++) {
                        // Check that the cell is in the grid!
                        int tmp_test = grid->test_cell(prim_id + delta);
                        if (tmp_test < 0)
                            continue;
                        Cell sec = grid->c[tmp_test];

                        // Define primary position
                        Float3 ppos = grid->p[j].pos;
                        #if PERIODIC
                        ppos -= grid->cell_sep(delta);
                        #endif

                        int end = sec.start + sec.np; // cell ID to stop the next loop at
                        if (j < end) end = j; // guarantee k < j to skip self-counts and count each pair only once, spare k < j check inside the loop

                        // This is the position of the particle as
                        // viewed from the secondary cell. Now loop over
                        // the particles in this secondary cell
                        for (int k = sec.start; k < end; k++) {
                            // end <= j thus k < j
                            // end <= sec.start + sec.np thus k is in the current cell
                            // Now we're considering these two particles!
                            Float3 dx = grid->p[k].pos - ppos;
                            Float norm2 = dx.norm2();
                            // Check if this is in the correct binning ranges
                            if ((norm2 < rmax_long2 && norm2 > rmin_long2) || (norm2 < rmax_cf2 && norm2 > rmin_cf2))
                                cnt2++;
                            else
                                continue;

                            // Now what do we want to do with the pair?
                            norm2 = sqrt(norm2);  // Now just radius
                            dx = dx / norm2; // normalize coordinate difference
                            // Find the long radial bin
                            int bin_long = floor((norm2 - rmin_long) / (rmax_long - rmin_long) * NBIN_LONG);
                            if ((bin_long >= 0) && (bin_long < NBIN_LONG)) { // if not out of bounds
                                // Accumulate the 4-pt correlation function
                                npcf[thread].add_4pcf(bin_long, pairs_i + j, pairs_i + k);
                                // Accumulate the mixed 3-pt correlation function
                                npcf[thread].add_3pcf_mixed(bin_long, pairs_i + j, grid->p[k].w);
                                npcf[thread].add_3pcf_mixed(bin_long, pairs_i + k, primary_w); // switch primary and secondary
                                #if (!PREVENT_TRIANGLES && !IGNORE_TRIANGLES)
                                // Exclude the mixed 3-pt correlation function self-counts
                                int bin_short = floor((norm2 - rmin_short) / (rmax_short - rmin_short) * NBIN_SHORT);
                                if ((bin_short >= 0) && (bin_short < NBIN_SHORT))
                                    npcf[thread].excl_3pcf_mixed(bin_long, bin_short, primary_w * grid->p[k].w * (primary_w + grid->p[k].w));
                                    // deal with jkj and kjk pseudo-triplets at once
                                #endif
                                // Accumulate long 2PCF
                                #if (PERIODIC)
                                if ((primary_w > 0) && (grid->p[k].w > 0)) // count only data-data pairs for 2PCF
                                #endif
                                npcf[thread].add_2pcf_long(bin_long, primary_w * grid->p[k].w);
                            }
                            // Find the fine 2pcf radial bin
                            int bin_cf = floor((norm2 - rmin_cf) / (rmax_cf - rmin_cf) * NBIN_CF);
                            // Get mu
                            #if (PERIODIC)
                            Float mu = dx.z; // cos(theta) with respect to z axis
                            #else
                            Float3 sumx = grid->p[k].pos + ppos;
                            Float sumx_norm = sumx.norm();
                            sumx = sumx / sumx_norm;
                            Float mu = sumx.dot(dx); // cos(theta) with respect to mean line of sight, to match gs4PCF
                            #endif
                            // Accumulate fine 2PCF
                            #if (PERIODIC)
                            if ((primary_w > 0) && (grid->p[k].w > 0)) // count only data-data pairs for 2PCF
                            #endif
                            if ((bin_cf >= 0) && (bin_cf < NBIN_CF)) // if not out of bounds
                                finepairs[thread].add(bin_cf, mu, grid->p[k].w * primary_w);
                        }  // Done with this secondary particle
                    }      // Done with this delta.z loop
                            // done with delta.y loop
            // done with delta.x loop
        }  // Done with this primary particle

    }  // Done with this primary cell (second loop), end of omp pragma

    #if (PERIODIC)
    // do analytical DR and RR counts for fine pairs
    for (int bin = 0; bin < NBIN_CF; bin++) {
        Float rmin_bin = rmin_cf + bin * (rmax_cf - rmin_cf) / NBIN_CF;
        Float rmax_bin = rmin_cf + (bin + 1) * (rmax_cf - rmin_cf) / NBIN_CF;
        Float w_prod = 0.5 * grid->sumw_neg * (2 * grid->sumw_pos + grid->sumw_neg) / pow(grid->max_boxsize, 3) * 4 * M_PI / 3 * (pow(rmax_bin, 3) - pow(rmin_bin, 3)) / MBIN_CF; // product of weights, contribution to a fine 2PCF bin; extra 1/2 because we count only k < j pairs
        for (int mubin = 0; mubin < MBIN_CF; mubin++) finepairs[0].add_raw(bin, mubin, w_prod); // add for each mu bin; to 0th thread because it always exists
    }
    // do analytical DR and RR counts for long pairs
    for (int bin = 0; bin < NBIN_LONG; bin++) {
        Float rmin_bin = rmin_long + bin * (rmax_long - rmin_long) / NBIN_LONG;
        Float rmax_bin = rmin_long + (bin + 1) * (rmax_long - rmin_long) / NBIN_LONG;
        Float w_prod = 0.5 * grid->sumw_neg * (2 * grid->sumw_pos + grid->sumw_neg) / pow(grid->max_boxsize, 3) * 4 * M_PI / 3 * (pow(rmax_bin, 3) - pow(rmin_bin, 3)); // product of weights, contribution to a long 2PCF bin; extra 1/2 because we count only k < j pairs
        npcf[0].add_2pcf_long(bin, w_prod); // add to 0th thread because it always exists
    }
    #endif

    powertime.Stop();

#ifndef OPENMP
#ifdef AVX
    printf(
        "\n# Time to compute required pairs (with AVX): %.2f\n\n",
        accpairs.Elapsed());
#else
    printf(
        "\n# Time to compute required pairs (no AVX): %.2f\n\n",
        accpairs.Elapsed());
#endif
#endif

    printf("# We counted  %lld pairs within [%f %f].\n", cnt, rmin_short, rmax_short);
    printf("# Average of %f pairs per primary particle.\n",
           (Float)cnt / grid->np);
    Float3 boxsize = grid->rect_boxsize;
    float expected = grid->np * (4 * M_PI / 3.0) *
                     (pow(rmax_short, 3.0) - pow(rmin_short, 3.0)) /
                     (boxsize.x * boxsize.y * boxsize.z);
    printf(
        "# We expected %1.0f pairs per primary particle, off by a factor of "
        "%f.\n",
        expected, cnt / (expected * grid->np));

    printf("# We counted  %lld triplets within [%f %f].\n", cnt3, rmin_short, rmax_short);

    printf("# We counted  %lld pairs within [%f %f].\n", cnt2, rmin_long, rmax_long);
    printf("# Average of %f pairs per primary particle.\n",
           (Float)cnt2 / grid->np);
    expected = grid->np * (4 * M_PI / 3.0) * (pow(rmax_long, 3.0) - pow(rmin_long, 3.0)) / (boxsize.x * boxsize.y * boxsize.z) / 2;
    printf(
        "# We expected %1.0f pairs per primary particle, off by a factor of "
        "%f.\n",
        expected, cnt2 / (expected * grid->np));

    // Detailed timing breakdown
    printf("\n# Accumulate Pairs: %6.3f s\n", accpairs.Elapsed());
    printf("# Compute Power: %6.3f s\n\n", powertime.Elapsed());

    delete[] pairs_i;

    return;
}

#endif
