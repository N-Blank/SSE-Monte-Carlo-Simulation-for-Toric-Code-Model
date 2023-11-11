#include "sampling.h"

/* 
 * function: sample 
 *  samples the current state in the simulation
 * 
 *  parameters:
 *      (int) n: bin number
 *      (int) t_idx: temperature index
 *      (toricode_system *) system: system to sample from
 *      (sse_state *) state: simulation state
 *      (sampled_quantities *) samples: store the samples
 */
void sample(int n, int t_idx, toricode_system *system, sse_state *state, sampled_quantities *samples, pcg32_random_t* rng) 
{
    double m1 = 0.0;
    double m2 = 0.0;
    double m4 = 0.0;
    double ms_tmp = 0.0;
    double ms = 0.0;
    double m2s = 0.0;
    double m4s = 0.0;
    double corr[system->L];
    double si[system->L];
    memset(corr, 0.0, system->L * sizeof(double));
    memset(si, 0.0, system->L * sizeof(double));

    // sample the first state 
    for (int i = 0; i < system->N; i++) {
        m1 += system->spin[i] * 0.5;
        ms_tmp += pow(- 1.0, i) * system->spin[i] * 0.5;
        // corr[i] += system->spin[0] * system->spin[i] * 0.25;
        // si[i] += system->spin[i] * 0.5;
    }
    m2 += pow(m1, 2.0);
    m4 += pow(m1, 4.0);
    ms += ms_tmp;
    m2s += pow(ms_tmp, 2.0);
    m4s += pow(ms_tmp, 4.0);

    // propagate the system to sample the staggered magnetization
    for (int p = 0; p < state->n; p++) {
        int b = (state->red_op_string[p] / 17) - 1;
        int a = state->red_op_string[p] % 17;
        if (a == 1) {
            ms_tmp += pow(-1.0, system->bond[b][0]) + pow(-1.0, system->bond[b][1]) + pow(-1.0, system->bond[b][2]) + pow(-1.0, system->bond[b][3] + pow(-1.0, system->bond[b][4]) + pow(-1.0, system->bond[b][5]));
        } 
        else if (a == 2) {
            ms_tmp += pow(-1.0, system->bond[b][0]) + pow(-1.0, system->bond[b][1]) + pow(-1.0, system->bond[b][2]) - pow(-1.0, system->bond[b][3] + pow(-1.0, system->bond[b][4]) + pow(-1.0, system->bond[b][5]));
        }
        else if (a == 3) {
            ms_tmp += pow(-1.0, system->bond[b][0]) + pow(-1.0, system->bond[b][1]) - pow(-1.0, system->bond[b][2]) + pow(-1.0, system->bond[b][3] + pow(-1.0, system->bond[b][4]) + pow(-1.0, system->bond[b][5]));
        }
        else if (a == 4) {
            ms_tmp += pow(-1.0, system->bond[b][0]) + pow(-1.0, system->bond[b][1]) - pow(-1.0, system->bond[b][2]) - pow(-1.0, system->bond[b][3] + pow(-1.0, system->bond[b][4]) + pow(-1.0, system->bond[b][5]));
        }
        else if (a == 5) {
            ms_tmp += pow(-1.0, system->bond[b][0]) - pow(-1.0, system->bond[b][1]) + pow(-1.0, system->bond[b][2]) + pow(-1.0, system->bond[b][3] + pow(-1.0, system->bond[b][4]) + pow(-1.0, system->bond[b][5]));
        }
        else if (a == 6) {
            ms_tmp += pow(-1.0, system->bond[b][0]) - pow(-1.0, system->bond[b][1]) + pow(-1.0, system->bond[b][2]) - pow(-1.0, system->bond[b][3] + pow(-1.0, system->bond[b][4]) + pow(-1.0, system->bond[b][5]));
        }
        else if (a == 7) {
            ms_tmp += pow(-1.0, system->bond[b][0]) - pow(-1.0, system->bond[b][1]) - pow(-1.0, system->bond[b][2]) + pow(-1.0, system->bond[b][3] + pow(-1.0, system->bond[b][4]) + pow(-1.0, system->bond[b][5]));
        }
        else if (a == 8) {
            ms_tmp += pow(-1.0, system->bond[b][0]) - pow(-1.0, system->bond[b][1]) - pow(-1.0, system->bond[b][2]) - pow(-1.0, system->bond[b][3] + pow(-1.0, system->bond[b][4]) + pow(-1.0, system->bond[b][5]));
        }
        else if (a == 9) {
            ms_tmp += -pow(-1.0, system->bond[b][0]) + pow(-1.0, system->bond[b][1]) + pow(-1.0, system->bond[b][2]) + pow(-1.0, system->bond[b][3] + pow(-1.0, system->bond[b][4]) + pow(-1.0, system->bond[b][5]));
        }
        else if (a == 10) {
            ms_tmp += -pow(-1.0, system->bond[b][0]) + pow(-1.0, system->bond[b][1]) + pow(-1.0, system->bond[b][2]) - pow(-1.0, system->bond[b][3] + pow(-1.0, system->bond[b][4]) + pow(-1.0, system->bond[b][5]));
        }
        else if (a == 11) {
            ms_tmp += -pow(-1.0, system->bond[b][0]) + pow(-1.0, system->bond[b][1]) - pow(-1.0, system->bond[b][2]) + pow(-1.0, system->bond[b][3] + pow(-1.0, system->bond[b][4]) + pow(-1.0, system->bond[b][5]));
        }
        else if (a == 12) {
            ms_tmp += -pow(-1.0, system->bond[b][0]) + pow(-1.0, system->bond[b][1]) - pow(-1.0, system->bond[b][2]) - pow(-1.0, system->bond[b][3] + pow(-1.0, system->bond[b][4]) + pow(-1.0, system->bond[b][5]));
        }
        else if (a == 13) {
            ms_tmp += -pow(-1.0, system->bond[b][0]) - pow(-1.0, system->bond[b][1]) + pow(-1.0, system->bond[b][2]) + pow(-1.0, system->bond[b][3] + pow(-1.0, system->bond[b][4]) + pow(-1.0, system->bond[b][5]));
        }
        else if (a == 14) {
            ms_tmp += -pow(-1.0, system->bond[b][0]) - pow(-1.0, system->bond[b][1]) + pow(-1.0, system->bond[b][2]) - pow(-1.0, system->bond[b][3] + pow(-1.0, system->bond[b][4]) + pow(-1.0, system->bond[b][5]));
        }
        else if (a == 15) {
            ms_tmp += -pow(-1.0, system->bond[b][0]) - pow(-1.0, system->bond[b][1]) - pow(-1.0, system->bond[b][2]) + pow(-1.0, system->bond[b][3] + pow(-1.0, system->bond[b][4]) + pow(-1.0, system->bond[b][5]));
        }
        else if (a == 16) {
            ms_tmp += -pow(-1.0, system->bond[b][0]) - pow(-1.0, system->bond[b][1]) - pow(-1.0, system->bond[b][2]) - pow(-1.0, system->bond[b][3] + pow(-1.0, system->bond[b][4]) + pow(-1.0, system->bond[b][5]));
        }
        ms += ms_tmp;
        m2s += pow(ms_tmp, 2.0);
        m4s += pow(ms_tmp, 4.0);

        // for (int i = 0; i < system->L; i++) {
        //     corr[i] += system->spin[0] * system->spin[i] * 0.25;
        //     si[i] += system->spin[i] * 0.5;
        // }
    }

    int norm = state->n + 1;
    // sample to the struct
    samples->n_bins[t_idx][n] += state->n;
    samples->n2_bins[t_idx][n] += state->n * state->n; 
    samples->m_bins[t_idx][n] += m1 / system->N;
    samples->m2_bins[t_idx][n] += m2 / pow(system->N, 2.0);
    samples->m4_bins[t_idx][n] += m4 / pow(system->N, 4.0);
    samples->ms_bins[t_idx][n] += ms / (system->N * norm);
    samples->m2s_bins[t_idx][n] += m2s / (pow(system->N, 2.0) * norm);
    samples->m4s_bins[t_idx][n] += m4s / (pow(system->N, 4.0) * norm);
    // for (int i = 0; i < system->L; i++) {
    //     samples->corr_bins[t_idx][n][i] += corr[i];
    //     samples->S_bins[t_idx][n] += pow(- 1.0, i) * corr[i] / system->L;
    // }

    // sample the kinetic coefficients

#ifdef KINETIC
    int *Sa = (int *) malloc((state->n + 1) * sizeof(int));
    int *Sb = (int *) malloc((state->n + 1) * sizeof(int));
    int *Ha = (int *) malloc(state->n * sizeof(int));
    int *Hb = (int *) malloc(state->n * sizeof(int));
    int Cab = 0;

#if defined(L_SS) || defined(L_SH)
    Sa[0] = 0.0;
    for (int a = samples->x + 1; a < system->L; a++) {
        Sa[0] += system->spin[a];
    }
#endif
#if defined(L_SS) || defined(L_HS)
    Sb[0] = 0.0;
    for (int b = samples->y + 1; b < system->L; b++) {
        Sb[0] += system->spin[b];
    }
#endif

    for (int p = 0; p < state->n; p++) {
        int bond = (state->red_op_string[p] / 17) - 1;
        int type = state->red_op_string[p] % 17;

        int change_a = 0;
        int change_b = 0;
        if (bond == samples->x && type != 0) {
            if (type == 1) {
                change_a = 8;
            } 
            else if (type == 2) {
                change_a = 4;
            }
            else if (type == 3) {
                change_a = 4;
            }
            else if (type == 4) {
                change_a = 0;
            }
            else if (type == 5) {
                change_a = 4;
            }
            else if (type == 6) {
                change_a = 0;
            }
            else if (type == 7) {
                change_a = 0;
            }
            else if (type == 8) {
                change_a = -4;
            }
            else if (type == 9) {
                change_a = 4;
            }
            else if (type == 10) {
                change_a = 0;
            }
            else if (type == 11) {
                change_a = 0;
            }
            else if (type == 12) {
                change_a = -4;
            }
            else if (type == 13) {
                change_a = 0;
            }
            else if (type == 14) {
                change_a = -4;
            }
            else if (type == 15) {
                change_a = -4;
            }
            else if (type == 16) {
                change_a = -8;
            }
        }
        if (bond == samples->y && type != 0) {
            if (type == 1) {
                change_b = 8;
            } 
            else if (type == 2) {
                change_b = 4;
            }
            else if (type == 3) {
                change_b = 4;
            }
            else if (type == 4) {
                change_b = 0;
            }
            else if (type == 5) {
                change_b = 4;
            }
            else if (type == 6) {
                change_b = 0;
            }
            else if (type == 7) {
                change_b = 0;
            }
            else if (type == 8) {
                change_b = -4;
            }
            else if (type == 9) {
                change_b = 4;
            }
            else if (type == 10) {
                change_b = 0;
            }
            else if (type == 11) {
                change_b = 0;
            }
            else if (type == 12) {
                change_b = -4;
            }
            else if (type == 13) {
                change_b = 0;
            }
            else if (type == 14) {
                change_b = -4;
            }
            else if (type == 15) {
                change_b = -4;
            }
            else if (type == 16) {
                change_b = -8;
            }
        }

#if defined(L_SS) || defined(L_SH)
        Sa[p + 1] = Sa[p] + change_a;
#endif
#if defined(L_SS) || defined(L_HS)
        Sb[p + 1] = Sb[p] + change_b;
#endif
#if defined(L_HH) || defined(L_HS)
        Ha[p] = (bond > samples->x) ? 1 : 0;
#endif
#if defined(L_HH) || defined(L_SH)
        Hb[p] = (bond > samples->y) ? 1 : 0;
#endif
#if defined(L_HH) 
        Cab += (bond > samples->y && bond > samples->x) ? 1 : 0;
#endif 
    }

    for (int samp = 0; samp < samples->max_samp; samp++) {
        // assing times to the operators
        double *tau = (double *) malloc((state->n + 2) * sizeof(double));
        tau[0] = 0.0;
        for (int i = 1; i <= state->n; i++) {
            tau[i] = pcg32_double_r(rng) * samples->beta_vals[t_idx];
        }
        tau[state->n + 1] = samples->beta_vals[t_idx];
        qsort(tau, state->n + 2, sizeof(double), compare);
        // for (int i = 1; i <= state->n; i++) {
        //     double tmp = pcg32_double_r(rng) * samples->beta_vals[t_idx];
        //     int j;
        //     for (j = i - 1; (j >= 1 && tau[j] > tmp); j--)
        //         tau[j + 1] = tau[j];
        //     tau[j + 1] = tmp;
        // }

        for (int k = 0; k < samples->k_max; k++) {
            double Ka[2] = {};
            double Kb[2] = {};
            double Ga[2] = {};
            double Gb[2] = {};

            for (int p = 0; p <= state->n; p++) {
#if defined(L_SS) || defined(L_SH)
                Ka[0] += (sin(samples->w_k[t_idx][k] * tau[p + 1]) - sin(samples->w_k[t_idx][k] * tau[p])) * 0.5 * Sa[p] / samples->w_k[t_idx][k];
                Ka[1] += (cos(samples->w_k[t_idx][k] * tau[p]) - cos(samples->w_k[t_idx][k] * tau[p + 1])) * 0.5 * Sa[p] / samples->w_k[t_idx][k];
#endif
#if defined(L_SS) || defined(L_HS)
                Kb[0] += (sin(samples->w_k[t_idx][k] * tau[p + 1]) - sin(samples->w_k[t_idx][k] * tau[p])) * 0.5 * Sb[p] / samples->w_k[t_idx][k];
                Kb[1] += (cos(samples->w_k[t_idx][k] * tau[p]) - cos(samples->w_k[t_idx][k] * tau[p + 1])) * 0.5 * Sb[p] / samples->w_k[t_idx][k];
#endif
#if defined(L_HH) || defined(L_HS) || defined(L_SH)
                if (p < state->n) {
#if defined(L_HH) || defined(L_HS)
                    Ga[0] += cos(samples->w_k[t_idx][k] * tau[p + 1]) * Ha[p];
                    Ga[1] += sin(samples->w_k[t_idx][k] * tau[p + 1]) * Ha[p];
#endif
#if defined(L_HH) || defined(L_SH)
                    Gb[0] += cos(samples->w_k[t_idx][k] * tau[p + 1]) * Hb[p];
                    Gb[1] += sin(samples->w_k[t_idx][k] * tau[p + 1]) * Hb[p];
#endif
                }
#endif
            }
            
#if defined(L_SS)
            samples->L_SS_bins[t_idx][n][k] += samples->w_k[t_idx][k] * (Ka[0] * Kb[0] + Ka[1] * Kb[1]) / (samples->beta_vals[t_idx] * samples->max_samp);
#endif 
#if defined(L_SH)
            if (state->n >= 1)
                samples->L_SH_bins[t_idx][n][k] += - samples->w_k[t_idx][k] * (Ka[0] * Gb[0] + Ka[1] * Gb[1]) / (samples->beta_vals[t_idx] * samples->max_samp);
#endif 
#if defined(L_HS)
            if (state->n >= 1)
                samples->L_HS_bins[t_idx][n][k] += - samples->w_k[t_idx][k] * (Ga[0] * Kb[0] + Ga[1] * Kb[1]) / (samples->beta_vals[t_idx] * samples->max_samp);
#endif
#if defined(L_HH)
            if (state->n >= 2)
                samples->L_HH_bins[t_idx][n][k] += samples->w_k[t_idx][k] * (Ga[0] * Gb[0] + Ga[1] * Gb[1] - Cab) / (samples->beta_vals[t_idx] * samples->max_samp);
#endif
        }
        free(tau);
    }
    free(Sa);
    free(Sb);
    free(Ha);
    free(Hb);
#endif // KINETIC
}

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
void normalize(long mc_cycles, sampled_quantities *samples, int N, double S) 
{
    for (int t_idx = 0; t_idx < samples->betas; t_idx++) {
        for (int n = 0; n < samples->bins; n++) {
            samples->n_bins[t_idx][n] /= mc_cycles;
            samples->n2_bins[t_idx][n] /= mc_cycles;
            samples->E_bins[t_idx][n] = - samples->n_bins[t_idx][n] 
                / (samples->beta_vals[t_idx] * N); /* remove the added constant to the Hamiltonian */
            
            samples->m_bins[t_idx][n] /= mc_cycles;
            samples->m2_bins[t_idx][n] /= mc_cycles;
            samples->m4_bins[t_idx][n] /= mc_cycles;
            samples->ms_bins[t_idx][n] /= mc_cycles;
            samples->m2s_bins[t_idx][n] /= mc_cycles;
            samples->m4s_bins[t_idx][n] /= mc_cycles;
            samples->m_sus_bins[t_idx][n] = samples->beta_vals[t_idx] 
                * N * (samples->m2_bins[t_idx][n] - samples->m_bins[t_idx][n] 
                * samples->m_bins[t_idx][n]);

            for (int i = 0; i < N; i++) {
                samples->corr_bins[t_idx][n][i] /= mc_cycles;
            }
            samples->S_bins[t_idx][n] /= mc_cycles;

            
            for (int k = 0; k < samples->k_max; k++) {
                samples->L_SS_bins[t_idx][n][k] /= mc_cycles;
                samples->L_HH_bins[t_idx][n][k] /= mc_cycles;
                samples->L_SH_bins[t_idx][n][k] /= mc_cycles;
                samples->L_HS_bins[t_idx][n][k] /= mc_cycles;
            }
            

            samples->n_mean[t_idx] += samples->n_bins[t_idx][n];
            samples->n2_mean[t_idx] += samples->n2_bins[t_idx][n];
            samples->E_mean[t_idx] += samples->E_bins[t_idx][n];
            samples->m_mean[t_idx] += samples->m_bins[t_idx][n];
            samples->m2_mean[t_idx] += samples->m2_bins[t_idx][n];
            samples->m4_mean[t_idx] += samples->m4_bins[t_idx][n];
            samples->ms_mean[t_idx] += samples->ms_bins[t_idx][n];
            samples->m2s_mean[t_idx] += samples->m2s_bins[t_idx][n];
            samples->m4s_mean[t_idx] += samples->m4s_bins[t_idx][n];
            samples->m_sus_mean[t_idx] += samples->m_sus_bins[t_idx][n];

            for (int i = 0; i < N; i++) {
                samples->corr_mean[t_idx][i] += samples->corr_bins[t_idx][n][i];
            }
            samples->S_mean[t_idx] += samples->S_bins[t_idx][n];

            
            for (int k = 0; k < samples->k_max; k++) {
                samples->L_SS_mean[t_idx][k] += samples->L_SS_bins[t_idx][n][k];
                samples->L_HH_mean[t_idx][k] += samples->L_HH_bins[t_idx][n][k];
                samples->L_SH_mean[t_idx][k] += samples->L_SH_bins[t_idx][n][k];
                samples->L_HS_mean[t_idx][k] += samples->L_HS_bins[t_idx][n][k];
            }
            
        }
        samples->n_mean[t_idx] /= samples->bins;
        samples->n2_mean[t_idx] /= samples->bins;
        samples->E_mean[t_idx] /= samples->bins;
        samples->m_mean[t_idx] /= samples->bins;
        samples->m2_mean[t_idx] /= samples->bins;
        samples->m4_mean[t_idx] /= samples->bins;
        samples->ms_mean[t_idx] /= samples->bins;
        samples->m2s_mean[t_idx] /= samples->bins;
        samples->m4s_mean[t_idx] /= samples->bins;
        samples->m_sus_mean[t_idx] /= samples->bins;

        for (int i = 0; i < N; i++) {
            samples->corr_mean[t_idx][i] /= samples->bins;
        }
        samples->S_mean[t_idx] /= samples->bins;

        
        for (int k = 0; k < samples->k_max; k++) {
            samples->L_SS_mean[t_idx][k] /= samples->bins;
            samples->L_HH_mean[t_idx][k] /= samples->bins;
            samples->L_SH_mean[t_idx][k] /= samples->bins;
            samples->L_HS_mean[t_idx][k] /= samples->bins;
        }
        

        for (int n = 0; n < samples->bins; n++) {
            samples->n_std[t_idx] += pow(samples->n_bins[t_idx][n] - samples->n_mean[t_idx], 2.0);
            samples->E_std[t_idx] += pow(samples->E_bins[t_idx][n] - samples->E_mean[t_idx], 2.0);
            samples->m_std[t_idx] += pow(samples->m_bins[t_idx][n] - samples->m_mean[t_idx], 2.0);
            samples->m2_std[t_idx] += pow(samples->m2_bins[t_idx][n] - samples->m2_mean[t_idx], 2.0);
            samples->m4_std[t_idx] += pow(samples->m4_bins[t_idx][n] - samples->m4_mean[t_idx], 2.0);
            samples->ms_std[t_idx] += pow(samples->ms_bins[t_idx][n] - samples->ms_mean[t_idx], 2.0);
            samples->m2s_std[t_idx] += pow(samples->m2s_bins[t_idx][n] - samples->m2s_mean[t_idx], 2.0);
            samples->m4s_std[t_idx] += pow(samples->m4s_bins[t_idx][n] - samples->m4s_mean[t_idx], 2.0);
            samples->m_sus_std[t_idx] += pow(samples->m_sus_bins[t_idx][n] - samples->m_sus_mean[t_idx], 2.0);

            for (int i = 0; i < N; i++) {
                samples->corr_std[t_idx][i] += pow(samples->corr_bins[t_idx][n][i] - samples->corr_mean[t_idx][i], 2.0);
            }
            samples->S_std[t_idx] += pow(samples->S_bins[t_idx][n] - samples->S_mean[t_idx], 2.0);

            
            for (int k = 0; k < samples->k_max; k++) {
                samples->L_SS_std[t_idx][k] += pow(samples->L_SS_bins[t_idx][n][k] - samples->L_SS_mean[t_idx][k], 2.0);
                samples->L_HH_std[t_idx][k] += pow(samples->L_HH_bins[t_idx][n][k] - samples->L_HH_mean[t_idx][k], 2.0);
                samples->L_SH_std[t_idx][k] += pow(samples->L_SH_bins[t_idx][n][k] - samples->L_SH_mean[t_idx][k], 2.0);
                samples->L_HS_std[t_idx][k] += pow(samples->L_HS_bins[t_idx][n][k] - samples->L_HS_mean[t_idx][k], 2.0);
            }
            
        }
        samples->n_std[t_idx] = sqrt(samples->n_std[t_idx] / samples->bins);
        samples->E_std[t_idx] = sqrt(samples->E_std[t_idx] / samples->bins);
        samples->m_std[t_idx] = sqrt(samples->m_std[t_idx] / samples->bins);
        samples->m2_std[t_idx] = sqrt(samples->m2_std[t_idx] / samples->bins);
        samples->m4_std[t_idx] = sqrt(samples->m4_std[t_idx] / samples->bins);
        samples->ms_std[t_idx] = sqrt(samples->ms_std[t_idx] / samples->bins);
        samples->m2s_std[t_idx] = sqrt(samples->m2s_std[t_idx] / samples->bins);
        samples->m4s_std[t_idx] = sqrt(samples->m4s_std[t_idx] / samples->bins);
        samples->m_sus_std[t_idx] = sqrt(samples->m_sus_std[t_idx] / samples->bins);

        for (int i = 0; i < N; i++) {
            samples->corr_std[t_idx][i] = sqrt(samples->corr_std[t_idx][i] / samples->bins);
        }
        samples->S_std[t_idx] = sqrt(samples->S_std[t_idx] / samples->bins);

        
        for (int k = 0; k < samples->k_max; k++) {
            samples->L_SS_std[t_idx][k] = sqrt(samples->L_SS_std[t_idx][k] / samples->bins);
            samples->L_HH_std[t_idx][k] = sqrt(samples->L_HH_std[t_idx][k] / samples->bins);
            samples->L_SH_std[t_idx][k] = sqrt(samples->L_SH_std[t_idx][k] / samples->bins);
            samples->L_HS_std[t_idx][k] = sqrt(samples->L_HS_std[t_idx][k] / samples->bins);
        }
        
    }
}

/*
 * function: compare
 * compares two numbers for the imaginary time assignments
 */
int compare( const void* num1, const void* num2)
{
    double a = *(double*) num1;  
    double b = *(double*) num2;  

    if(a > b)
    {  
        return 1;  
    }  
    else if(a < b)  
    {  
        return -1;  
    }  
    return 0;  
}
