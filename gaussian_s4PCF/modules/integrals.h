// Rewritten integrals.h code for grid_covariance.cpp (originally from Alex Wiegand) to parallelize and compute integrands to a given quad of particles
#include "parameters.h"
#include "correlation_function.h"
#include "cell_utilities.h"
#include "jackknife_weights.h"

#ifndef INTEGRALS_H
#define INTEGRALS_H

class Integrals{
private:
    CorrelationFunction *cf12, *cf13, *cf24;
    Float rmin,rmax,rmin_long,rmax_long,mumin,mumax,dmu; //Ranges in r and mu
    Float *r_high, *r_low, *r_high_long, *r_low_long; // Max and min of each radial bin
    int nbin, nbin_long, mbin;
    int Nbin, Nbinpairs, N4;
    Float *c4; // Array to accumulate integral
    Float *norm4; // Array to store normalizations for integral
    JK_weights *JK12, *JK23, *JK34; // RR counts and jackknife weights
    char* out_file;
    bool box,rad=0; // Flags to decide whether we have a periodic box + if we have a radial correlation function only
    int I1, I2, I3, I4; // indices for which fields to use for each particle

    uint64 *binct4; // Array to accumulate bin counts

public:
    Integrals(){};

    Integrals(Parameters *par, CorrelationFunction *_cf12, CorrelationFunction *_cf13, CorrelationFunction *_cf24, JK_weights *_JK12, JK_weights *_JK23, JK_weights *_JK34, int _I1, int _I2, int _I3, int _I4){
        cf12 = new CorrelationFunction(_cf12);
        cf13 = new CorrelationFunction(_cf13);
        cf24 = new CorrelationFunction(_cf24);
        JK12 = _JK12;
        JK23 = _JK23;
        JK34 = _JK34;
        I1=_I1;
        I2=_I2;
        I3=_I3;
        I4=_I4;
        init(par);
    }
    void init(Parameters *par){
        nbin = par->nbin; // number of radial bins
        nbin_long = par->nbin_long; // number of long radial bins
        mbin = par->mbin; // number of mu bins
        Nbin = nbin*mbin; // number of short side bins
        Nbinpairs = Nbin*(Nbin+1)/2; // number of unique short side bin pairs
        N4 = nbin_long*Nbinpairs; // total bumber of bins in squeezed 4PCF
        out_file = par->out_file; // output directory

        int ec=0;
        // Initialize the binning
        ec+=posix_memalign((void **) &c4, PAGE, sizeof(double)*N4);
        ec+=posix_memalign((void **) &norm4, PAGE, sizeof(double)*N4);
        ec+=posix_memalign((void **) &binct4, PAGE, sizeof(uint64)*N4);

        assert(ec==0);
        reset();

        box=par->perbox;

        rmax=par->rmax;
        rmin=par->rmin;

        rmax_long=par->rmax_long;
        rmin_long=par->rmin_long;

        mumax=par->mumax;
        mumin=par->mumin;

        r_high = par->radial_bins_high;
        r_low = par->radial_bins_low;

        r_high_long = par->radial_bins_high_long;
        r_low_long = par->radial_bins_low_long;

        dmu=(mumax-mumin)/mbin;

        rad=mbin==1&&dmu==1.;
    }

    ~Integrals() {
        free(c4);
        free(binct4);
    }

    void reset(){
        for (int j=0; j<N4; j++) {
            c4[j]=0;
            binct4[j] = 0;
        }
    }

    inline int getbin(Float r, Float mu){
        // Linearizes 2D indices
        // First define which r bin we are in;
        int which_bin = -1; // default if outside bins
        for(int i=0;i<nbin;i++){
            if((r>r_low[i])&&(r<r_high[i])){
                which_bin=i;
                break;
            }
            if((i==nbin-1)&&(r>r_high[i])){
                which_bin=nbin; // if above top bin
            }
        }
        return which_bin*mbin + floor((mu-mumin)/dmu);
    }

    inline int getbin_pair(int bin1, int bin2){
        int i = bin1, j = bin2;
        if (j < i) {
            i = bin2;
            j = bin1;
        }
        return j + Nbin*i - i * (i-1) / 2;
    }

    inline int getbin_long(Float r, Float mu){
        // Linearizes 2D indices
        // First define which r bin we are in;
        int which_bin = -1; // default if outside bins
        for(int i=0;i<nbin;i++){
            if((r>r_low_long[i])&&(r<r_high_long[i])){
                which_bin=i;
                break;
            }
            if((i==nbin-1)&&(r>r_high[i])){
                which_bin=nbin; // if above top bin
            }
        }
        return which_bin;
    }

    inline void fourth(const Particle* pi_list, const int* prim_ids, const int pln, const Particle pj, const Particle pk, const Particle pl, const int pj_id, const int pk_id, const int pl_id, const double prob){
        // Accumulates the four point integral C4.
        // First define variables
        Particle pi;
        Float rjl_mag, rjl_mu, rkl_mag, rkl_mu, rij_mag, rij_mu, rik_mag, rik_mu, ril_mag, ril_mu, rjk_mag, rjk_mu, c4v, tmp_weight;
        int bin_jl, tmp_full_bin;
        cleanup_l(pl.pos,pk.pos,rkl_mag,rkl_mu);
        Float xi_kl = cf24->xi(rkl_mag, rkl_mu); // should be cf34 in general but that does not exist yet
        cleanup_l(pl.pos,pj.pos,rjl_mag,rjl_mu);
        bin_jl = getbin(rjl_mag,rjl_mu);
        if ((bin_jl<0)||(bin_jl>=Nbin)) return; // if not in correct bin
        Float xi_jl = cf24->xi(rjl_mag, rjl_mu); // j-l correlation
        cleanup_l(pk.pos,pj.pos,rjk_mag,rjk_mu);
        Float xi_jk = cf24->xi(rjk_mag, rjk_mu); // should be cf23 in general but that does not exist yet

        for(int i=0;i<pln;i++){ // Iterate over particle in pi_list
            if(((prim_ids[i]==pj_id)&&(I1==I2))||((prim_ids[i]==pk_id)&&(I1==I3))) continue; // skip ij/ik self counts
            if(((prim_ids[i]==pl_id)&&(I1==I4))||((pj_id==pl_id)&&(I2==I4))||((pk_id==pl_id)&&(I3==I4))) continue; // don't self-count

            pi = pi_list[i];
            Float wijk = pi.w*pj.w*pk.w;
            tmp_weight = wijk*pl.w; // product of weights, w_i*w_j*w_k*w_l

            cleanup_l(pi.pos, pj.pos, rij_mag, rij_mu);
            int bin_ij = getbin_long(rij_mag, rij_mu);
            if ((bin_ij<0)||(bin_ij>=nbin_long)) continue; // if not in correct bin
            Float xi_ij = cf12->xi(rij_mag, rij_mu);
            cleanup_l(pi.pos, pk.pos, rik_mag, rik_mu);
            int bin_ik = getbin(rik_mag, rik_mu);
            if ((bin_ik<0)||(bin_ik>=Nbin)) continue; // if not in correct bin
            Float xi_ik = cf13->xi(rik_mag, rik_mu);
            cleanup_l(pi.pos, pl.pos, ril_mag, ril_mu);
            Float xi_il = cf12->xi(ril_mag, ril_mu); // should be cf14 in general, but that does not exist yet
            tmp_full_bin = bin_ij*Nbinpairs + getbin_pair(bin_ik, bin_jl);

            // Now compute the integral;
            c4v = tmp_weight/prob*(xi_ik*xi_jl + xi_ij*xi_kl + xi_il*xi_jk); // gaussian 4PCF
            // Add to local counts
            c4[tmp_full_bin]+=c4v;
            norm4[tmp_full_bin]+=tmp_weight/prob;
            binct4[tmp_full_bin]++;
        }
    }

public:
    void cleanup_l(Float3 p1,Float3 p2,Float& norm,Float& mu){
        Float3 pos=p1-p2;
        norm = pos.norm();
        // If the input correlation function had only radial information and no mu bins
        // fill the dummy mu bins by 0.5
        if(rad){
            mu=0.5;
        }
        else{
#ifndef PERIODIC
            Float3 los=p1+p2; // No 1/2 as normalized anyway below
            mu = fabs(pos.dot(los)/norm/los.norm());
#else
            // In the periodic case use z-direction for mu
            mu = fabs(pos.z/norm);
#endif
        }
    }
public:
    void sum_ints(Integrals* ints) {
        // Add the values accumulated in ints to the corresponding internal sums
        for(int i=0; i<N4; i++) {
            c4[i]+=ints->c4[i];
            norm4[i]+=ints->norm4[i];
            binct4[i]+=ints->binct4[i];
        }
    }

    void frobenius_difference_sum(Integrals* ints, int n_loop, Float &frobC4){
        // Add the values accumulated in ints to the corresponding internal sums and compute the Frobenius norm difference between integrals
        Float n_loops = (Float)n_loop;
        Float self_c4=0, diff_c4=0;
        // Compute Frobenius norms and sum integrals
        for(int i=0; i<N4; i++) {
            self_c4+=pow(c4[i]/n_loops,2.);
            diff_c4+=pow(c4[i]/n_loops-(c4[i]+ints->c4[i])/(n_loops+1.),2.);
        }
        diff_c4=sqrt(diff_c4);
        self_c4=sqrt(self_c4);
        // Return percent difference
        frobC4=100.*(diff_c4/self_c4);
        }

    void sum_total_counts(uint64& acc4){
        // Add local counts to total bin counts in acc4
        for(int i=0; i<N4; i++) {
            acc4+=binct4[i];
        }
    }

    void normalize(){
        // Normalize the accumulated integrals
        for(int i=0; i<N4; i++) {
            c4[i]/=norm4[i];
        }
    }
    void save_counts(uint64 quad_counts){
        // Print the counts for each integral (used for combining the estimates outside of C++)
        // This is the number of counts used in each loop [always the same]
        char counts_file[1000];
        snprintf(counts_file, sizeof counts_file, "%sCovMatricesAll/total_counts_n%d_m%d_%d%d,%d%d.txt",out_file,nbin,mbin,I1,I2,I3,I4);
        FILE * CountsFile = fopen(counts_file,"w");
        fprintf(CountsFile,"%llu\n",quad_counts);

        fflush(NULL);
        fclose(CountsFile);
    }


    void save_integrals(char* suffix, bool save_all) {
    /* Print integral outputs to file.
        * In txt files {c2,c3,c4,RR}_n{nbin}_m{mbin}.txt there are lists of the outputs of c2,c3,c4 and RR_a that are already normalized and multiplied by combinatoric factors. The n and m strings specify the number of n and m bins present.
        */
        // Create output files

        char c4name[1000];
        snprintf(c4name, sizeof c4name, "%sCovMatricesAll/c4_n%d_m%d_%d%d,%d%d_%s.txt", out_file, nbin, mbin, I1,I2,I3,I4,suffix);

        FILE * C4File = fopen(c4name,"w"); // for c4 part of integral

        for(int i=0; i<N4; i++) {
            fprintf(C4File,"%le\t",c4[i]);
        }
        fprintf(C4File,"\n");

        fflush(NULL);

        // Close open files
        fclose(C4File);

        if(save_all==1){
            char binname[1000];
            snprintf(binname,sizeof binname, "%sCovMatricesAll/binct_c4_n%d_m%d_%d%d,%d%d_%s.txt",out_file, nbin,mbin,I1,I2,I3,I4,suffix);
            FILE * BinFile = fopen(binname,"w");

            for(int i=0; i<N4; i++) {
                fprintf(BinFile,"%llu\t",binct4[i]);
            }
            fprintf(BinFile,"\n");
            fclose(BinFile);
        }
    }

};

#endif
