// jackknife_weights.h - this contains c++ functions to read-in pair counts from files.

#ifndef JACKKNIFE_WEIGHTS_H
#define JACKKNIFE_WEIGHTS_H

// This class stores the RR count weights for a given jackknife
class JK_weights{
public:
    Float* RR_pair_counts; // houses the weighted pair counts summed over jackknife regions.
    int nbins; // total number of bins
    
public: 
    
    void copy(JK_weights *JK){
        // Copy JK_weights object
        nbins=JK->nbins;
        RR_pair_counts = (Float *)malloc(sizeof(Float)*nbins);
        for(int i=0;i<nbins;i++){
            RR_pair_counts[i]=JK->RR_pair_counts[i];
        }
    }
    
    void rescale(Float norm1, Float norm2){
        // Rescale the RR pair counts by a factor (N_gal1/N_rand1)*(N_gal2/N_rand2)
        Float rescale_factor = norm1*norm2;
        printf("Rescaling RR pair counts by a factor (N_gal_1/N_rand_1)*(N_gal2/N_rand2) = %.1e\n",1./rescale_factor);
        for(int i=0;i<nbins;i++){
            RR_pair_counts[i]/=rescale_factor;
        }
    }
    
public:
    ~JK_weights() {
        // The destructor
        free(RR_pair_counts);
        return;
    }
    
    JK_weights(){};
    
    // Assignment operator creation
    JK_weights& operator=(const JK_weights& jk);
    
    JK_weights(Parameters *par, int index1, int index2){
        
        // This reads in weights for each jackknife region for each bin from file.
        // File should have bins in space-separated columns and bins in rows, with the jackknife number in the first column.
        // NB: The indexing is defined as INDEX = JACKKNIFE_ID * NBINS + BIN_ID
        // index1 and index2 define which random set of particles to use here
        // If Jackknife directive is not set, we only read in RR pair counts here
        
        nbins = par->nbin*par->mbin; // define number of bins in total
        char line[1000000];
        FILE *fp2;
        char *RR_file;
        
        if((index1==1)&&(index2==1)) RR_file = par->RR_bin_file;
        else if((index1==2)&&(index2==2)) RR_file = par->RR_bin_file2;
        else RR_file = par->RR_bin_file12;
        fp2 = fopen(RR_file,"r");
        
        if (fp2==NULL){
            fprintf(stderr,"RR bin count file %s not found\n",RR_file);
            abort();
        }
        fprintf(stderr,"\nReading RR bin count file '%s'\n",RR_file);
        
        int ec2=0;
        ec2+=posix_memalign((void **) &RR_pair_counts, PAGE, sizeof(Float)*nbins);
        assert(ec2==0);
        
        int index=0;
        while (fgets(line,5000,fp2)!=NULL){
            // Select required lines in file
            if (line[0]=='#') continue;
            if (line[0]=='\n') continue;
            RR_pair_counts[index]=atof(line);
            index++;
        }
        assert(index==nbins);
        printf("Read in RR pair counts successfully.\n");
    }
};
  
#endif
