#ifndef SAMPLING_H
#define SAMPLING_H

#include "../sse/sse.h"
#include <time.h>

/*
 * struct: sampled_quantities
 *  contains all of the sampled quantities during the simulation
 *  
 *  quantities:
 *      energy
 *      specific heat
 *      magnetization
 *      staggered magnetization
 *      magnetic uniform susceptibility
 *      binder parameter
 */
typedef struct sampled_quantities 
{
    int L;
    int bins;
    int betas;
    double *beta_vals;

    double **n_bins;
    double **n2_bins;
    double **E_bins;
    double **m_bins;
    double **m2_bins;
    double **m4_bins;
    double **ms_bins;
    double **m2s_bins;
    double **m4s_bins;
    double **m_sus_bins;

    double *n_mean;
    double *n2_mean;
    double *n_std;

    double *E_mean;
    double *E_std;

    double *m_mean;
    double *m_std;
    double *m2_mean;
    double *m2_std;
    double *m4_mean;
    double *m4_std;
    double *ms_mean;
    double *ms_std;
    double *m2s_mean;
    double *m2s_std;
    double *m4s_mean;
    double *m4s_std;
    double *m_sus_mean;
    double *m_sus_std;

    double ***corr_bins;
    double **corr_mean;
    double **corr_std;

    double **S_bins;
    double *S_mean;
    double *S_std;


    int k_max;
    int x;
    int y;
    int max_samp;
    double **w_k;
    double ***L_SS_bins;
    double **L_SS_mean;
    double **L_SS_std;
    double ***L_HH_bins;
    double **L_HH_mean;
    double **L_HH_std;
    double ***L_SH_bins;
    double **L_SH_mean;
    double **L_SH_std;
    double ***L_HS_bins;
    double **L_HS_mean;
    double **L_HS_std;
} sampled_quantities;

/* 
 * function: init_samples
 *  initializes and allocates memory for the samples quantities struct
 * 
 *  parameters:
 *      (double *) beta_vals: array of the temperatures
 *      (int) len_beta: number of temperatures
 *      (int) n_bins: number of bins
 *      (int) d: number of dimensions of the system
 *      (int) L: number of unit lattices in the system
 *      (sampled_quantities *) samples: struct to inilialize 
 *      (int) max_samp: number of samples for the conductance
 */
void init_samples(double *beta_vals, int len_beta, int n_bins, int L,
    sampled_quantities *samples, int max_samp);

/* 
 * function: sample 
 *  samples the currect state in the simulation
 * 
 *  parameters:
 *      (int) n: bin number
 *      (int) t_idx: temperature index
 *      (toricode_system *) system: system to sample from
 *      (sse_state *) state: simulation state
 *      (sampled_quantities *) samples: store the samples
 */
void sample(int n, int t_idx, toricode_system *system, sse_state *state, 
    sampled_quantities *samples, pcg32_random_t* rng);

/* 
 * function: normalize 
 *  normalize the samples during the binning phase
 * 
 *  parameters:
 *      (long) mc_cycles: MCS for sampling
 *      (sampled_quantities *) samples: sampled quantities
 *      (int) N: number of particles
 *      (double) S: quantum spin number
 */
void normalize(long mc_cycles, sampled_quantities *samples, int N, double S);

/*
 * function: free_samples
 *  frees the allocated memory in the sampled_quantities struct
 * 
 *  parameters:
 *      (sampled_quantities *) samples: struct to free
 */
void free_samples(sampled_quantities *samples);

/*
 * function: compare
 * compares two numbers for the imaginary time assignments
 */
int compare( const void* num1, const void* num2);

#endif // SAMPLING_H
