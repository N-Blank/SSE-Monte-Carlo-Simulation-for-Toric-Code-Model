#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>
#include <mpi.h>



#include "sse/sse.h"
#include "sampling/sampling.h"
#include "io/io.h"

#ifndef TESTING
#define SEED (u_int64_t) time(NULL)
#else
#define SEED (u_int64_t) 2
#endif

// global variables for the system and simulation
int L;
double S;

long therm_cycles;
long mc_cycles;
int n_slots;
int n_bins;
int n_threads;

double *beta_vals;
int len_beta;

sampled_quantities *samples[100];

// variables for debugging
double *cpu_time;
int *n_locals;

void simulate(int myid, int t_id, char *vtx_file)
{
    toricode_system *system = (toricode_system *) malloc(sizeof(toricode_system));
    sse_state *state = (sse_state *) malloc(sizeof(sse_state));

    init_toricode_system(L, S, system);
    init_sse_state(SEED * t_id, system, state);

    pcg32_random_t rng;
    pcg32_srandom_r(&rng, (SEED * t_id) ^ (intptr_t)&rng, (SEED * t_id));

    read_vtx_info(vtx_file, &(state->vtx_type), &(state->n_diagrams));

    for (int t_idx = 0; t_idx < len_beta; t_idx++) {
        clock_t start_clock = clock();
        double beta = beta_vals[t_idx];
        if (t_idx > 0) { reset_sse_state(system, state); }
        for (long t = 0; t < therm_cycles; t++) {
            diag_update(beta, system, state, &rng);

            create_vtx_list(system, state);
            loffdiagonal_update(system, state, &rng);

            ajust_cutoff(state, t);
        }
        
        #pragma omp parallel for
        for (int n = 0; n < n_bins; n++) {
            #pragma omp parallel for
            for (long t = 0; t < mc_cycles; t++) {
                diag_update(beta, system, state, &rng);
                
                create_vtx_list(system, state);
                loffdiagonal_update(system, state, &rng);

                sample(n, t_idx, system, state, samples[myid], &rng);
            }
        }

        clock_t end_clock = clock();

        time_t t = time(NULL);
        char *buff = ctime(&t);
        buff[strcspn(buff, "\n")] = 0;

        cpu_time[t_id - 1] = ((double) (end_clock - start_clock)) / CLOCKS_PER_SEC;
        n_locals[t_id - 1] = state->n_locals;

        {
            double max = cpu_time[0];
            int max_id = 0;
            int avg_n_locals = 0;
            for (int i = 0; i < n_threads; i++) {
                if (max < cpu_time[i]) {
                    max = cpu_time[i];
                    max_id = i;
                }
                avg_n_locals += n_locals[i];
            }
            avg_n_locals /= n_threads;

            printf("%s | slotsID: %d | beta: %.4lf | n_locals: %d | max_time: %.5lfs (T_id: %d) \n", 
                buff, myid, beta, avg_n_locals, max / n_threads, max_id);
            fflush(stdout);
        }
    }
    free_memory(system, state);
    free(state);
    free(system);
}

int main(int argc, char **argv)
{
    if (argc < 3) {
        printf("Please provide the output file names for the program to work. \n");
        printf("Usage: mpi -np nslots %s n_threads vtx_name.txt \n", argv[0]);
        exit(1);
    }

    read_inputs(&L, &S, &therm_cycles, &mc_cycles, &n_bins, &beta_vals, &len_beta);
    n_threads = atoi(argv[1]);
    n_slots = 10;

    int myid, numprocs;
    sampled_quantities *global_samples;
    global_samples = (sampled_quantities *) malloc(sizeof(sampled_quantities));
    init_samples(beta_vals, len_beta, n_bins, L, global_samples, 1);
    int size = global_samples->betas;
    int k_max = global_samples->k_max;
    int sizeL = size * L;
    int sizek_max = size * k_max;

    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD,&myid);

    samples[myid] = (sampled_quantities *) malloc(sizeof(sampled_quantities));
    int max_samp = 1;
    if (argc == 7) {
        max_samp = atoi(argv[6]);
    }

    init_samples(beta_vals, len_beta, n_bins, L, samples[myid], max_samp);
    
    cpu_time = (double *) malloc(sizeof(double) * n_threads);
    n_locals = (int *) malloc(sizeof(int) * n_threads);

    time_t t = time(NULL);
    printf(" -- Starting SSE simulation %d of the spin-S toricode model -- \n", myid);
    printf(" L: %d | S: %.1lf", L, S);
    printf("n_slots: %d | slots_id: %d |  n_threads: %d | therm_cycles: %ld | mc_cycles: %ld | n_bins: %d \n", 
            numprocs, myid, n_threads, therm_cycles, mc_cycles, n_bins);

    printf("   Simulation %d started at: %s ", myid, ctime(&t));
    printf("\n");
    fflush(stdout);
    clock_t start_clock = clock();
    
    int t_id = 0;
    simulate(myid, t_id + 1, argv[2]);


    printf("\n");
    normalize(mc_cycles, samples[myid], 2 * pow(L, 2), L, S);
   
    MPI_Reduce(&samples[myid]->n_std[0], &global_samples->n_std[0], size, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&samples[myid]->E_std[0], &global_samples->E_std[0], size, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&samples[myid]->Ov_std[0], &global_samples->Ov_std[0], size, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&samples[myid]->m_std[0], &global_samples->m_std[0], size, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&samples[myid]->m2_std[0], &global_samples->m2_std[0], size, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&samples[myid]->m4_std[0], &global_samples->m4_std[0], size, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&samples[myid]->ms_std[0], &global_samples->ms_std[0], size, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&samples[myid]->m2s_std[0], &global_samples->m2s_std[0], size, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&samples[myid]->m4s_std[0], &global_samples->m4s_std[0], size, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&samples[myid]->m_sus_std[0], &global_samples->m_sus_std[0], size, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&samples[myid]->S_std[0], &global_samples->S_std[0], size, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&samples[myid]->n_mean[0], &global_samples->n_mean[0], size, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&samples[myid]->n2_mean[0], &global_samples->n2_mean[0], size, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&samples[myid]->E_mean[0], &global_samples->E_mean[0], size, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&samples[myid]->Ov_mean[0], &global_samples->Ov_mean[0], size, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&samples[myid]->m_mean[0], &global_samples->m_mean[0], size, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&samples[myid]->m2_mean[0], &global_samples->m2_mean[0], size, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&samples[myid]->m4_mean[0], &global_samples->m4_mean[0], size, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&samples[myid]->ms_mean[0], &global_samples->ms_mean[0], size, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&samples[myid]->m2s_mean[0], &global_samples->m2s_mean[0], size, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&samples[myid]->m4s_mean[0], &global_samples->m4s_mean[0], size, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&samples[myid]->m_sus_mean[0], &global_samples->m_sus_mean[0], size, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&samples[myid]->S_mean[0], &global_samples->S_mean[0], size, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    
    if (myid == 0){
        for (int t_idx = 0; t_idx < samples[0]->betas; t_idx++){
            global_samples->n_std[t_idx] = sqrt(global_samples->n_std[t_idx] / n_slots);
            global_samples->E_std[t_idx] = sqrt(global_samples->E_std[t_idx] / n_slots);
            global_samples->m_std[t_idx] = sqrt(global_samples->m_std[t_idx] / n_slots);
            global_samples->m2_std[t_idx] = sqrt(global_samples->m2_std[t_idx] / n_slots);
            global_samples->m4_std[t_idx] = sqrt(global_samples->m4_std[t_idx] / n_slots);
            global_samples->ms_std[t_idx] = sqrt(global_samples->ms_std[t_idx] / n_slots);
            global_samples->m2s_std[t_idx] = sqrt(global_samples->m2s_std[t_idx] / n_slots);
            global_samples->m4s_std[t_idx] = sqrt(global_samples->m4s_std[t_idx] / n_slots);
            global_samples->m_sus_std[t_idx] = sqrt(global_samples->m_sus_std[t_idx] / n_slots);
            global_samples->S_std[t_idx] = sqrt(global_samples->S_std[t_idx] / n_slots);
            global_samples->n_mean[t_idx] = global_samples->n_mean[t_idx] / n_slots;
            global_samples->E_mean[t_idx] = global_samples->E_mean[t_idx] / n_slots;
            global_samples->m_mean[t_idx] = global_samples->m_mean[t_idx] / n_slots;
            global_samples->m2_mean[t_idx] = global_samples->m2_mean[t_idx] / n_slots;
            global_samples->m4_mean[t_idx] = global_samples->m4_mean[t_idx] / n_slots;
            global_samples->ms_mean[t_idx] = global_samples->ms_mean[t_idx] / n_slots;
            global_samples->m2s_mean[t_idx] = global_samples->m2s_mean[t_idx] / n_slots;
            global_samples->m4s_mean[t_idx] = global_samples->m4s_mean[t_idx] / n_slots;
            global_samples->m_sus_mean[t_idx] = global_samples->m_sus_mean[t_idx] / n_slots;
            global_samples->S_mean[t_idx] = global_samples->S_mean[t_idx] / n_slots;

            for (int i = 0; i < L; i++) {
                global_samples->corr_mean[t_idx][i] = global_samples->corr_mean[t_idx][i] / n_slots;
                global_samples->corr_std[t_idx][i] = sqrt(global_samples->corr_std[t_idx][i] / n_slots);
            }

            for (int k = 0; k < samples[0]->k_max; k++) {
                global_samples->L_SS_mean[t_idx][k] = global_samples->L_SS_mean[t_idx][k] / n_slots;
                global_samples->L_HH_mean[t_idx][k] = global_samples->L_HH_mean[t_idx][k] / n_slots;
                global_samples->L_SH_mean[t_idx][k] = global_samples->L_SH_mean[t_idx][k] / n_slots;
                global_samples->L_HS_mean[t_idx][k] = global_samples->L_HS_mean[t_idx][k] / n_slots;
                global_samples->L_SS_std[t_idx][k] = sqrt(global_samples->L_SS_std[t_idx][k] / n_slots);
                global_samples->L_HH_std[t_idx][k] = sqrt(global_samples->L_HH_std[t_idx][k] / n_slots);
                global_samples->L_SH_std[t_idx][k] = sqrt(global_samples->L_SH_std[t_idx][k] / n_slots);
                global_samples->L_HS_std[t_idx][k] = sqrt(global_samples->L_HS_std[t_idx][k] / n_slots);
            }
        }
        clock_t end_clock = clock();
        double cpu_time_used = ((double) (end_clock - start_clock)) / (CLOCKS_PER_SEC * n_threads);
        free(cpu_time);
        free(n_locals);
        
        printf("Simulation %d finished in %.5lfs \n", myid, cpu_time_used);
        printf(" -- Writing simulation results to file -- \n");
        fflush(stdout);
        char *file_name = write_outputs(global_samples, L, S, therm_cycles, mc_cycles, cpu_time_used, n_threads);
        printf(" -- Results written with success to file: %s -- \n", file_name);
        fflush(stdout);
        free(file_name);
    }
    free_samples(samples[myid]);
    MPI_Finalize();
}