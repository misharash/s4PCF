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
    int nbin, mbin;
    Float rmin,rmax,mumin,mumax,dmu; //Ranges in r and mu
    Float *r_high, *r_low; // Max and min of each radial bin
    Float *c4; // Array to accumulate integral
    JK_weights *JK12, *JK23, *JK34; // RR counts and jackknife weights
    char* out_file;
    bool box,rad=0; // Flags to decide whether we have a periodic box + if we have a radial correlation function only
    int I1, I2, I3, I4; // indices for which fields to use for each particle

    uint64 *binct, *binct3, *binct4; // Arrays to accumulate bin counts

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
        mbin = par->mbin; // number of mu bins
        out_file = par->out_file; // output directory

        int ec=0;
        // Initialize the binning
        ec+=posix_memalign((void **) &c4, PAGE, sizeof(double)*nbin*mbin*nbin*mbin);
        ec+=posix_memalign((void **) &binct, PAGE, sizeof(uint64)*nbin*mbin);
        ec+=posix_memalign((void **) &binct3, PAGE, sizeof(uint64)*nbin*mbin*nbin*mbin);
        ec+=posix_memalign((void **) &binct4, PAGE, sizeof(uint64)*nbin*mbin*nbin*mbin);

        assert(ec==0);
        reset();

        box=par->perbox;

        rmax=par->rmax;
        rmin=par->rmin;

        mumax=par->mumax;
        mumin=par->mumin;

        r_high = par->radial_bins_high;
        r_low = par->radial_bins_low;

        dmu=(mumax-mumin)/mbin;

        rad=mbin==1&&dmu==1.;
    }

    ~Integrals() {
        free(c4);
        free(binct);
        free(binct3);
        free(binct4);
    }

    void reset(){
        for (int j=0; j<nbin*mbin; j++) {
            binct[j] = 0;
        }
        for (int j=0; j<nbin*mbin*nbin*mbin; j++) {
            c4[j]=0;
            binct3[j] = 0;
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

    inline void fourth(const Particle* pi_list, const int* prim_ids, const int pln, const Particle pj, const Particle pk, const Particle pl, const int pj_id, const int pk_id, const int pl_id, const int* bin_ij, const Float* wijk, const Float* xi_ik, const double prob){
        // Accumulates the four point integral C4.
        // First define variables
        Particle pi;
        Float rjl_mag, rjl_mu, rkl_mag, rkl_mu, c4v, xi_jl, tmp_weight;
        int tmp_bin, tmp_full_bin, max_bin = nbin*mbin;
        cleanup_l(pl.pos,pk.pos,rkl_mag,rkl_mu);
        
        tmp_bin = getbin(rkl_mag,rkl_mu);

        if ((tmp_bin<0)||(tmp_bin>=max_bin)) return; // if not in correct bin
        cleanup_l(pl.pos,pj.pos,rjl_mag,rjl_mu);
        xi_jl = cf24->xi(rjl_mag, rjl_mu); // j-l correlation

        for(int i=0;i<pln;i++){ // Iterate over particle in pi_list
            if(wijk[i]==-1) continue; // skip incorrect bins / ij self counts
            if(((prim_ids[i]==pl_id)&&(I1==I4))||((pj_id==pl_id)&&(I2==I4))||((pk_id==pl_id)&&(I3==I4))) continue; // don't self-count

            pi = pi_list[i];
            tmp_weight = wijk[i]*pl.w; // product of weights, w_i*w_j*w_k*w_l

            // Now compute the integral;
            c4v = tmp_weight/prob*2.*xi_ik[i]*xi_jl; // with xi_ik*xi_jl = xi_il*xi_jk symmetry factor
            // Compute jackknife weight tensor:
            tmp_full_bin = bin_ij[i]*mbin*nbin+tmp_bin;
            // Add to local counts
            c4[tmp_full_bin]+=c4v;
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
        for(int i=0;i<nbin*mbin;i++){
            binct[i]+=ints->binct[i];
        }
        for(int i=0;i<nbin*mbin;i++){
            for(int j=0;j<nbin*mbin;j++){
                c4[i*nbin*mbin+j]+=ints->c4[i*nbin*mbin+j];
                binct3[i*nbin*mbin+j]+=ints->binct3[i*nbin*mbin+j];
                binct4[i*nbin*mbin+j]+=ints->binct4[i*nbin*mbin+j];
            }
        }
    }

    void frobenius_difference_sum(Integrals* ints, int n_loop, Float &frobC2, Float &frobC3, Float &frobC4){
        // Add the values accumulated in ints to the corresponding internal sums and compute the Frobenius norm difference between integrals
        Float n_loops = (Float)n_loop;
        Float self_c2=0, diff_c2=0;
        Float self_c3=0, diff_c3=0;
        Float self_c4=0, diff_c4=0;
        // Compute Frobenius norms and sum integrals
        for(int i=0;i<nbin*mbin;i++){
            for(int j=0;j<nbin*mbin;j++){
                self_c4+=pow(c4[i*nbin*mbin+j]/n_loops,2.);
                diff_c4+=pow(c4[i*nbin*mbin+j]/n_loops-(c4[i*nbin*mbin+j]+ints->c4[i*nbin*mbin+j])/(n_loops+1.),2.);
                }
            }
        self_c2=sqrt(self_c2);
        diff_c2=sqrt(diff_c2);
        diff_c3=sqrt(diff_c3);
        diff_c4=sqrt(diff_c4);
        self_c3=sqrt(self_c3);
        self_c4=sqrt(self_c4);
        // Return percent difference
        frobC2=100.*(diff_c2/self_c2);
        frobC3=100.*(diff_c3/self_c3);
        frobC4=100.*(diff_c4/self_c4);
        }

    void sum_total_counts(uint64& acc2, uint64& acc3, uint64& acc4){
        // Add local counts to total bin counts in acc2-4
        for (int i=0; i<nbin*mbin; i++) {
            acc2+=binct[i];
            for (int j=0; j<nbin*mbin; j++) {
                acc3+=binct3[i*nbin*mbin+j];
                acc4+=binct4[i*nbin*mbin+j];
            }
        }
    }
    void normalize(Float norm1, Float norm2, Float norm3, Float norm4, Float n_pairs, Float n_triples, Float n_quads){
        // Normalize the accumulated integrals (partly done by the normalising probabilities used from the selected cubes)
        // n_pair etc. are the number of PARTICLE pairs etc. attempted (not including rejected cells, but including pairs which don't fall in correct bin ranges)
        // To avoid recomputation
        double corrf2 = norm1*norm2; // correction factor for densities of random points
        double corrf3 = corrf2*norm3;
        double corrf4 = corrf3*norm4;
        for(int i = 0; i<nbin*mbin;i++){
            for(int j=0;j<nbin*mbin;j++){
                c4[i*nbin*mbin+j]/=(n_quads*corrf4);
            }
        }

        // Further normalize by RR counts from corrfunc
        for(int i=0; i<nbin*mbin;i++){
            Float Ra_i = JK12->RR_pair_counts[i];
            for(int j=0;j<nbin*mbin;j++){
                Float Rab4=Ra_i*JK34->RR_pair_counts[j];
                c4[i*nbin*mbin+j]/=Rab4;
            }
        }
    }
    void save_counts(uint64 pair_counts,uint64 triple_counts,uint64 quad_counts){
        // Print the counts for each integral (used for combining the estimates outside of C++)
        // This is the number of counts used in each loop [always the same]
        char counts_file[1000];
        snprintf(counts_file, sizeof counts_file, "%sCovMatricesAll/total_counts_n%d_m%d_%d%d,%d%d.txt",out_file,nbin,mbin,I1,I2,I3,I4);
        FILE * CountsFile = fopen(counts_file,"w");
        fprintf(CountsFile,"%llu\n",pair_counts);
        fprintf(CountsFile,"%llu\n",triple_counts);
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

        for(int i=0;i<nbin*mbin;i++){
            for(int j=0;j<nbin*mbin;j++){
                fprintf(C4File,"%le\t",c4[i*nbin*mbin+j]);
            }
            fprintf(C4File,"\n");
        }

        fflush(NULL);

        // Close open files
        fclose(C4File);

        if(save_all==1){
            char binname[1000];
            snprintf(binname,sizeof binname, "%sCovMatricesAll/binct_c4_n%d_m%d_%d%d,%d%d_%s.txt",out_file, nbin,mbin,I1,I2,I3,I4,suffix);
            FILE * BinFile = fopen(binname,"w");

            char bin3name[1000];
            snprintf(bin3name,sizeof bin3name, "%sCovMatricesAll/binct_c3_n%d_m%d_%d,%d%d_%s.txt",out_file, nbin,mbin,I2,I1,I3,suffix);
            FILE * Bin3File = fopen(bin3name,"w");

            char bin2name[1000];
            snprintf(bin2name,sizeof bin2name, "%sCovMatricesAll/binct_c2_n%d_m%d_%d%d_%s.txt",out_file, nbin,mbin,I1,I2,suffix);
            FILE * Bin2File = fopen(bin2name,"w");

            for (int j=0;j<nbin*mbin;j++){
                fprintf(Bin2File,"%llu\n",binct[j]);
            }

            for(int i=0;i<nbin*mbin;i++){
                for(int j=0;j<nbin*mbin;j++){
                    fprintf(BinFile,"%llu\t",binct4[i*nbin*mbin+j]);
                    fprintf(Bin3File,"%llu\t",binct3[i*nbin*mbin+j]);
                    }
                fprintf(BinFile,"\n");
                fprintf(Bin3File,"\n");
            }
            fclose(BinFile);
            fclose(Bin2File);
            fclose(Bin3File);
        }
    }

};

#endif
