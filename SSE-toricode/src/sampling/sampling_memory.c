#include "sampling.h"

#define BUFFER_SIZE 256

/* 
 * function: init_samples
 *  initializes and allocates memory for the samples quantities struct
 * 
 *  parameters:
 *      (double *) beta_vals: array of the temperatures
 *      (int) len_beta: number of temperatures
 *      (int) n_bins: number of bins
 *      (int) L: number of unit lattices in the system
 *      (sampled_quantities *) samples: struct to inilialize 
 *      (int) max_samp: number of samples for the conductance
 */
void init_samples(double *beta_vals, int len_beta, int n_bins, int L, struct sampled_quantities *samples, int max_samp) 
{
    samples->bins = n_bins;
    samples->betas = len_beta;
    samples->beta_vals = beta_vals;
    samples->L = L;
    samples->k_max = 0;
    samples->x = 0;
    samples->y = 0;
    samples->max_samp = max_samp;
    

    samples->n_bins = (double **) malloc(len_beta * sizeof(double *));
    samples->n2_bins = (double **) malloc(len_beta * sizeof(double *));
    samples->E_bins = (double **) malloc(len_beta * sizeof(double *));
    samples->m_bins = (double **) malloc(len_beta * sizeof(double *));
    samples->m2_bins = (double **) malloc(len_beta * sizeof(double *));
    samples->m4_bins = (double **) malloc(len_beta * sizeof(double *));
    samples->ms_bins = (double **) malloc(len_beta * sizeof(double *));
    samples->m2s_bins = (double **) malloc(len_beta * sizeof(double *));
    samples->m4s_bins = (double **) malloc(len_beta * sizeof(double *));
    samples->m_sus_bins = (double **) malloc(len_beta * sizeof(double *));
    for (int i = 0; i < len_beta; i++) {
        samples->n_bins[i] = (double *) malloc(n_bins * sizeof(double));
        samples->n2_bins[i] = (double *) malloc(n_bins * sizeof(double));
        samples->E_bins[i] = (double *) malloc(n_bins * sizeof(double));
        samples->m_bins[i] = (double *) malloc(n_bins * sizeof(double));
        samples->m2_bins[i] = (double *) malloc(n_bins * sizeof(double));
        samples->m4_bins[i] = (double *) malloc(n_bins * sizeof(double));
        samples->ms_bins[i] = (double *) malloc(n_bins * sizeof(double));
        samples->m2s_bins[i] = (double *) malloc(n_bins * sizeof(double));
        samples->m4s_bins[i] = (double *) malloc(n_bins * sizeof(double));
        samples->m_sus_bins[i] = (double *) malloc(n_bins * sizeof(double));

        memset(samples->n_bins[i], 0.0, n_bins * sizeof(double));
        memset(samples->n2_bins[i], 0.0, n_bins * sizeof(double));
        memset(samples->E_bins[i], 0.0, n_bins * sizeof(double));
        memset(samples->m_bins[i], 0.0, n_bins * sizeof(double));
        memset(samples->m2_bins[i], 0.0, n_bins * sizeof(double));
        memset(samples->m4_bins[i], 0.0, n_bins * sizeof(double));
        memset(samples->ms_bins[i], 0.0, n_bins * sizeof(double));
        memset(samples->m2s_bins[i], 0.0, n_bins * sizeof(double));
        memset(samples->m4s_bins[i], 0.0, n_bins * sizeof(double));
        memset(samples->m_sus_bins[i], 0.0, n_bins * sizeof(double));
    }

    samples->n_mean = (double *) malloc(len_beta * sizeof(double));
    samples->n_std = (double *) malloc(len_beta * sizeof(double));
    samples->n2_mean = (double *) malloc(len_beta * sizeof(double));

    samples->E_mean = (double *) malloc(len_beta * sizeof(double));
    samples->E_std = (double *) malloc(len_beta * sizeof(double));

    samples->m_mean = (double *) malloc(len_beta * sizeof(double));
    samples->m_std = (double *) malloc(len_beta * sizeof(double));
    samples->m2_mean = (double *) malloc(len_beta * sizeof(double));
    samples->m2_std = (double *) malloc(len_beta * sizeof(double));
    samples->m4_mean = (double *) malloc(len_beta * sizeof(double));
    samples->m4_std = (double *) malloc(len_beta * sizeof(double));
    samples->ms_mean = (double *) malloc(len_beta * sizeof(double));
    samples->ms_std = (double *) malloc(len_beta * sizeof(double));
    samples->m2s_mean = (double *) malloc(len_beta * sizeof(double));
    samples->m2s_std = (double *) malloc(len_beta * sizeof(double));
    samples->m4s_mean = (double *) malloc(len_beta * sizeof(double));
    samples->m4s_std = (double *) malloc(len_beta * sizeof(double));
    samples->m_sus_mean = (double *) malloc(len_beta * sizeof(double));
    samples->m_sus_std = (double *) malloc(len_beta * sizeof(double));

    memset(samples->n_mean, 0.0, len_beta * sizeof(double));
    memset(samples->n_std, 0.0, len_beta * sizeof(double));
    memset(samples->n2_mean, 0.0, len_beta * sizeof(double));
    memset(samples->E_mean, 0.0, len_beta * sizeof(double));
    memset(samples->E_std, 0.0, len_beta * sizeof(double));
    memset(samples->m_mean, 0.0, len_beta * sizeof(double));
    memset(samples->m_std, 0.0, len_beta * sizeof(double));
    memset(samples->m2_mean, 0.0, len_beta * sizeof(double));
    memset(samples->m2_std, 0.0, len_beta * sizeof(double));
    memset(samples->m4_mean, 0.0, len_beta * sizeof(double));
    memset(samples->m4_std, 0.0, len_beta * sizeof(double));
    memset(samples->ms_mean, 0.0, len_beta * sizeof(double));
    memset(samples->ms_std, 0.0, len_beta * sizeof(double));
    memset(samples->m2s_mean, 0.0, len_beta * sizeof(double));
    memset(samples->m2s_std, 0.0, len_beta * sizeof(double));
    memset(samples->m4s_mean, 0.0, len_beta * sizeof(double));
    memset(samples->m4s_std, 0.0, len_beta * sizeof(double));
    memset(samples->m_sus_mean, 0.0, len_beta * sizeof(double));
    memset(samples->m_sus_std, 0.0, len_beta * sizeof(double));

    samples->corr_bins = (double ***) malloc(len_beta * sizeof(double **));
    samples->corr_mean = (double **) malloc(len_beta * sizeof(double *));
    samples->corr_std = (double **) malloc(len_beta * sizeof(double *));
    samples->S_bins = (double **) malloc(len_beta * sizeof(double *));
    samples->S_mean = (double *) malloc(len_beta * sizeof(double));
    samples->S_std = (double *) malloc(len_beta * sizeof(double));
    for (int i = 0; i < len_beta; i++) {
        samples->corr_bins[i] = (double **) malloc(n_bins * sizeof(double *));
        samples->corr_mean[i] = (double *) malloc(L * sizeof(double));
        samples->corr_std[i] = (double *) malloc(L * sizeof(double));

        samples->S_bins[i] = (double *) malloc(n_bins * sizeof(double));

        for (int j = 0; j < n_bins; j++) {
            samples->corr_bins[i][j] = (double *) malloc(L * sizeof(double));

            memset(samples->corr_bins[i][j], 0.0, L * sizeof(double));
        }
        
        memset(samples->corr_mean[i], 0.0, L * sizeof(double));
        memset(samples->corr_std[i], 0.0, L * sizeof(double));

        memset(samples->S_bins[i], 0.0, n_bins * sizeof(double));

    }

    memset(samples->S_mean, 0.0, len_beta * sizeof(double));
    memset(samples->S_std, 0.0, len_beta * sizeof(double));

#ifdef KINETIC
    char buffer[BUFFER_SIZE];
    FILE *fp;
    fp = fopen("matsubara.in", "r");
    if (fp != NULL) {
        fgets(buffer, BUFFER_SIZE, fp);
        sscanf(buffer, "%d ", &(samples->k_max));

        fgets(buffer, BUFFER_SIZE, fp);
        sscanf(buffer, "%d, %d", &(samples->x), &(samples->y));
    } else {
        printf("Error opening the matsubara.in file. Check if the file exists. \n");
        printf("Setting the maximum Matsubara frequency to 5. \n");
        samples->k_max = 5;
        printf("Setting x and y for measuring the perturbation to the middle of the chain. \n");
        samples->x = L / 2 - 1;
        samples->y = samples->x;
    }

    samples->L_SS_bins = (double ***) malloc(len_beta * sizeof(double **));
    samples->L_SS_mean = (double **) malloc(len_beta * sizeof(double *));
    samples->L_SS_std = (double **) malloc(len_beta * sizeof(double *));
    samples->L_HH_bins = (double ***) malloc(len_beta * sizeof(double **));
    samples->L_HH_mean = (double **) malloc(len_beta * sizeof(double *));
    samples->L_HH_std = (double **) malloc(len_beta * sizeof(double *));
    samples->L_SH_bins = (double ***) malloc(len_beta * sizeof(double **));
    samples->L_SH_mean = (double **) malloc(len_beta * sizeof(double *));
    samples->L_SH_std = (double **) malloc(len_beta * sizeof(double *));
    samples->L_HS_bins = (double ***) malloc(len_beta * sizeof(double **));
    samples->L_HS_mean = (double **) malloc(len_beta * sizeof(double *));
    samples->L_HS_std = (double **) malloc(len_beta * sizeof(double *));
    samples->w_k = (double **) malloc(len_beta * sizeof(double *));
    for (int i = 0; i < len_beta; i++) {
        samples->L_SS_bins[i] = (double **) malloc(n_bins * sizeof(double *));
        samples->L_SS_mean[i] = (double *) malloc(samples->k_max * sizeof(double));
        samples->L_SS_std[i] = (double *) malloc(samples->k_max * sizeof(double));
        samples->L_HH_bins[i] = (double **) malloc(n_bins * sizeof(double *));
        samples->L_HH_mean[i] = (double *) malloc(samples->k_max * sizeof(double));
        samples->L_HH_std[i] = (double *) malloc(samples->k_max * sizeof(double));
        samples->L_SH_bins[i] = (double **) malloc(n_bins * sizeof(double *));
        samples->L_SH_mean[i] = (double *) malloc(samples->k_max * sizeof(double));
        samples->L_SH_std[i] = (double *) malloc(samples->k_max * sizeof(double));
        samples->L_HS_bins[i] = (double **) malloc(n_bins * sizeof(double *));
        samples->L_HS_mean[i] = (double *) malloc(samples->k_max * sizeof(double));
        samples->L_HS_std[i] = (double *) malloc(samples->k_max * sizeof(double));
        samples->w_k[i] = (double *) malloc(samples->k_max * sizeof(double));

        for (int j = 0; j < n_bins; j++) {
            samples->L_SS_bins[i][j] = (double *) malloc(samples->k_max * sizeof(double));
            samples->L_HH_bins[i][j] = (double *) malloc(samples->k_max * sizeof(double));
            samples->L_SH_bins[i][j] = (double *) malloc(samples->k_max * sizeof(double));
            samples->L_HS_bins[i][j] = (double *) malloc(samples->k_max * sizeof(double));

            memset(samples->L_SS_bins[i][j], 0.0, samples->k_max * sizeof(double));
            memset(samples->L_HH_bins[i][j], 0.0, samples->k_max * sizeof(double));
            memset(samples->L_SH_bins[i][j], 0.0, samples->k_max * sizeof(double));
            memset(samples->L_HS_bins[i][j], 0.0, samples->k_max * sizeof(double));
        }

        memset(samples->L_SS_mean[i], 0.0, samples->k_max * sizeof(double));
        memset(samples->L_SS_std[i], 0.0, samples->k_max * sizeof(double));
        memset(samples->L_HH_mean[i], 0.0, samples->k_max * sizeof(double));
        memset(samples->L_HH_std[i], 0.0, samples->k_max * sizeof(double));
        memset(samples->L_SH_mean[i], 0.0, samples->k_max * sizeof(double));
        memset(samples->L_SH_mean[i], 0.0, samples->k_max * sizeof(double));
        memset(samples->L_HS_std[i], 0.0, samples->k_max * sizeof(double));
        memset(samples->L_HS_std[i], 0.0, samples->k_max * sizeof(double));
        memset(samples->w_k[i], 0.0, samples->k_max * sizeof(double));
    }

    for (int i = 0; i < len_beta; i++) {
        for (int k = 1; k < samples->k_max + 1; k++) {
            samples->w_k[i][k - 1] = 2 * M_PI * k / samples->beta_vals[i];
        }
    }
#endif // KINETIC
}

/*
 * function: free_samples
 *  frees the allocated memory in the sampled_quantities struct
 * 
 *  parameters:
 *      (sampled_quantities *) samples: struct to free
 */
void free_samples(sampled_quantities *samples) 
{
    for (int i = 0; i < samples->betas; i++) {
        free(samples->n_bins[i]);
        free(samples->n2_bins[i]);
        free(samples->E_bins[i]);
        free(samples->m_bins[i]);
        free(samples->m2_bins[i]);
        free(samples->m4_bins[i]);
        free(samples->ms_bins[i]);
        free(samples->m2s_bins[i]);
        free(samples->m4s_bins[i]);
        free(samples->m_sus_bins[i]);
    }
    free(samples->n_bins);
    free(samples->n2_bins);
    free(samples->E_bins);
    free(samples->m_bins);
    free(samples->m2_bins);
    free(samples->m4_bins);
    free(samples->ms_bins);
    free(samples->m2s_bins);
    free(samples->m4s_bins);
    free(samples->m_sus_bins);
    free(samples->n_mean);
    free(samples->n_std);
    free(samples->n2_mean);
    free(samples->E_mean);
    free(samples->E_std);
    free(samples->m_mean);
    free(samples->m_std);
    free(samples->m2_mean);
    free(samples->m2_std);
    free(samples->m4_mean);
    free(samples->m4_std);
    free(samples->ms_mean);
    free(samples->ms_std);
    free(samples->m2s_mean);
    free(samples->m2s_std);
    free(samples->m4s_mean);
    free(samples->m4s_std);
    free(samples->m_sus_mean);
    free(samples->m_sus_std);

    for (int i = 0; i < samples->betas; i++) {
        for (int j = 0; j < samples->bins; j++) {
            free(samples->corr_bins[i][j]);
        }
        free(samples->corr_bins[i]);
        free(samples->corr_mean[i]);
        free(samples->corr_std[i]);

        free(samples->S_bins[i]);
    }
    free(samples->corr_bins);
    free(samples->corr_mean);
    free(samples->corr_std);
    
    free(samples->S_bins);
    free(samples->S_mean);
    free(samples->S_std);

#ifdef KINETIC
    for (int i = 0; i < samples->betas; i++) {
        for (int j = 0; j < samples->bins; j++) {
            free(samples->L_SS_bins[i][j]);
            free(samples->L_HH_bins[i][j]);
            free(samples->L_SH_bins[i][j]);
            free(samples->L_HS_bins[i][j]);
        }   
        free(samples->L_SS_bins[i]);
        free(samples->L_SS_mean[i]);
        free(samples->L_SS_std[i]);
        free(samples->L_HH_bins[i]);
        free(samples->L_HH_mean[i]);
        free(samples->L_HH_std[i]);
        free(samples->L_SH_bins[i]);
        free(samples->L_SH_mean[i]);
        free(samples->L_SH_std[i]);
        free(samples->L_HS_bins[i]);
        free(samples->L_HS_mean[i]);
        free(samples->L_HS_std[i]);
        free(samples->w_k[i]);
    }
    free(samples->L_SS_bins);
    free(samples->L_SS_mean);
    free(samples->L_SS_std);
    free(samples->L_HH_bins);
    free(samples->L_HH_mean);
    free(samples->L_HH_std);
    free(samples->L_SH_bins);
    free(samples->L_SH_mean);
    free(samples->L_SH_std);
    free(samples->L_HS_bins);
    free(samples->L_HS_mean);
    free(samples->L_HS_std);
    free(samples->w_k);
#endif // KINETIC
}
