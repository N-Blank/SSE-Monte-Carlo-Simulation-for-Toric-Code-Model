#include "io.h"

#define BUFFER_SIZE 256

/*
 * function: read_inputs
 *  reads input file for the simulation and stores them in the 
 * parameters of the function
 *  
 *  parameters:
 *      (int *) L: length of the system
 *      (double *) S: spin quantum number
 *      (long *) therm_cycles: MCS for thermalization
 *      (long *) mc_cycles: MCS for sampling
 *      (int *) n_bins: number of sampling bins
 *      (double **) beta_vals: array to store simulation temperatures
 *      (int *) len_betas: number of temperatures 
 */
void read_inputs(int *L, double *S,  long *therm_cycles, long *mc_cycles,
 int *n_bins, double **beta_vals, int *len_beta) 
{
    char buffer[BUFFER_SIZE];
    FILE *input_file;

    input_file = fopen("read.in", "r");
    if (input_file != NULL) {
        fgets(buffer, BUFFER_SIZE, input_file);
        sscanf(buffer, "%d, %lf ", L, S);
        
        fgets(buffer, BUFFER_SIZE, input_file);
        sscanf(buffer, "%ld, %ld, %d ", therm_cycles, mc_cycles, n_bins);
    } else {
        printf("Error opening the read.i file. Check if the file exists. \n");
        exit(1);
    }
    fclose(input_file);

    input_file = fopen("beta.in", "r");
    if (input_file != NULL) {
        fgets(buffer, BUFFER_SIZE, input_file);
        sscanf(buffer, "%d ", len_beta);
        (*beta_vals) = (double *) malloc((*len_beta) * sizeof(double));
        
        for (int i = 0; i < (*len_beta); i++) {
            fgets(buffer, BUFFER_SIZE, input_file);
            sscanf(buffer, "%lf ", &((*beta_vals)[i]));
        }
    } else {
        printf("Error opening the beta.in file. Check if the file exists. \n");
        exit(1);
    }
    fclose(input_file);
}

/*
 * function: read_vtx_info
 *  reads vertex information for the simulation and 
 * stores it in the parameters of the function
 *  
 *  parameters:
 *      (char *) file_name: file name
 *      (vtx_element **) d: dimension
 *      (int *) n_diagrams: number of diagrams
 */
void read_vtx_info(char *file_name, vtx_element **vtx, int *n_diagrams) 
{
    FILE *vtx_file;
    vtx_file = fopen(file_name, "r");

    if (vtx_file != NULL) {
        fscanf(vtx_file, "%d \n", n_diagrams);
        (*vtx) = (vtx_element *) malloc((*n_diagrams) * sizeof(vtx_element));

        for (int i = 0; i < (*n_diagrams); i++) {
            fscanf(vtx_file, "%d \n", &((*vtx)[i].indx));
            fscanf(vtx_file, "%d \n", &((*vtx)[i].type));
            fscanf(vtx_file, "%lf \n", &((*vtx)[i].H));
            fscanf(vtx_file, "%d %d %d %d %d %d %d %d %d %d %d %d \n", &((*vtx)[i].spin[0]), &((*vtx)[i].spin[1]), &((*vtx)[i].spin[2]), &((*vtx)[i].spin[3]), &((*vtx)[i].spin[4]), &((*vtx)[i].spin[5]), &((*vtx)[i].spin[6]), &((*vtx)[i].spin[7]), &((*vtx)[i].spin[8]), &((*vtx)[i].spin[9]), &((*vtx)[i].spin[10]), &((*vtx)[i].spin[11]));
            fscanf(vtx_file, "%d %d \n", &((*vtx)[i].new_vtx_type[0]), &((*vtx)[i].new_vtx_type[1]));
        }
    } else {
        printf("Error opening the %s file. \n", file_name);
        exit(1);
    }

    fclose(vtx_file);
}

/*
 * function: write_outputs
 *  writes the simulation outputs to file
 *  
 *  parameters:
 *      (sampled_quantities *) samples: sampled quantities
 * during the simulation
 *  returns:
 *      (char *) file_name: name of save file
 */
char *write_outputs(sampled_quantities *samples, 
    int L, double S, long therm_cycles, long mc_cycles,
     double cpu_time_used, int n_threads) 
{
    char *file_name = (char *) malloc(BUFFER_SIZE * sizeof(char));
    sprintf(file_name, "L%d_toricode_S%g.csv", L, S);
    FILE *output_file;
    output_file = fopen(file_name, "w");
    if (output_file != NULL) {
        fprintf(output_file, "L,S\n");
        fprintf(output_file, "%d,%lf\n", L, S);
        
        fprintf(output_file, "therm_cycles,mc_cycles,n_bins\n");
        fprintf(output_file, "%ld,%ld,%d \n", therm_cycles, mc_cycles, samples->bins);
        
        fprintf(output_file, "cpu_time,n_threads\n");
        fprintf(output_file, "%lf,%d\n", cpu_time_used, n_threads);


#if defined(L_SS) && defined(L_HH) && defined(L_SH)
        char cond[5] = "full";
#else 
    #if defined(L_SS) && defined(L_HH)
        char cond[5] = "diag";
    #elif defined(L_SH)
        char cond[5] = "offd";
    #elif defined(L_HH) && defined(L_SH)
        char cond[5] = "hhsh"; 
    #elif defined(L_SS) && defined(L_SH)
        char cond[5] = "sssh";
    #elif defined(L_HH)
        char cond[5] = "hh";
    #elif defined(L_SS)
        char cond[5] = "ss";
    #else 
        char cond[5] = "."; 
    #endif
#endif
        fprintf(output_file, "n_betas,n_k,x,y,kinetic\n");
        fprintf(output_file, "%d,%d,%d,%d,%s\n", samples->betas, samples->k_max, samples->x, samples->y, cond);

        fprintf(output_file, "beta,n,n2,n_std,E,E_std,Ov,Ov_std,m,m_std,m2,m2_std,m4,m4_std,ms,ms_std,m2s,m2s_std,m4s,m4s_std,sus,sus_std,S_mean,S_std\n");
        for (int t_idx = 0; t_idx < samples->betas; t_idx++) {
            fprintf(output_file, "%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf\n", 
            samples->beta_vals[t_idx], 
            samples->n_mean[t_idx],
            samples->n2_mean[t_idx], 
            samples->n_std[t_idx],
            samples->E_mean[t_idx],
            samples->E_std[t_idx],
            samples->Ov_mean[t_idx],
            samples->Ov_std[t_idx],
            samples->m_mean[t_idx],
            samples->m_std[t_idx],
            samples->m2_mean[t_idx],
            samples->m2_std[t_idx],
            samples->m4_mean[t_idx],
            samples->m4_std[t_idx],
            samples->ms_mean[t_idx],
            samples->ms_std[t_idx],
            samples->m2s_mean[t_idx],
            samples->m2s_std[t_idx],
            samples->m4s_mean[t_idx],
            samples->m4s_std[t_idx],
            samples->m_sus_mean[t_idx],
            samples->m_sus_std[t_idx],
            samples->S_mean[t_idx],
            samples->S_std[t_idx]);
        }
        
        // for (int t_idx = 0; t_idx < samples->betas; t_idx++) {
        //     fprintf(output_file, "beta\n");
        //     fprintf(output_file, "%lf\n", samples->beta_vals[t_idx]);
        //     fprintf(output_file, "corr_mean,corr_std\n");
        //     for (int i = 0; i < L; i++) {
        //         fprintf(output_file, "%lf,%lf\n", samples->corr_mean[t_idx][i], samples->corr_std[t_idx][i]);
        //     }
        }


#ifdef L_SS
        for (int t_idx = 0; t_idx < samples->betas; t_idx++) {
            fprintf(output_file, "beta\n");
            fprintf(output_file, "%lf\n", samples->beta_vals[t_idx]);
            fprintf(output_file, "w_k,L_SS_mean,L_SS_std\n");
            for (int k = 0; k < samples->k_max; k++) {
                fprintf(output_file, "%lf,%lf,%lf\n", samples->w_k[t_idx][k], samples->L_SS_mean[t_idx][k], samples->L_SS_std[t_idx][k]);
            }
        }
#endif // L_SS
#ifdef L_HH
        for (int t_idx = 0; t_idx < samples->betas; t_idx++) {
            fprintf(output_file, "beta\n");
            fprintf(output_file, "%lf\n", samples->beta_vals[t_idx]);
            fprintf(output_file, "w_k,L_HH_mean,L_HH_std\n");
            for (int k = 0; k < samples->k_max; k++) {
                fprintf(output_file, "%lf,%lf,%lf\n", samples->w_k[t_idx][k], samples->L_HH_mean[t_idx][k], samples->L_HH_std[t_idx][k]);
            }
        }
#endif // L_HH
#ifdef L_SH
        for (int t_idx = 0; t_idx < samples->betas; t_idx++) {
            fprintf(output_file, "beta\n");
            fprintf(output_file, "%lf\n", samples->beta_vals[t_idx]);
            fprintf(output_file, "w_k,L_SH_mean,L_SH_std\n");
            for (int k = 0; k < samples->k_max; k++) {
                fprintf(output_file, "%lf,%lf,%lf\n", samples->w_k[t_idx][k], samples->L_SH_mean[t_idx][k], samples->L_SH_std[t_idx][k]);
            }
        }

        for (int t_idx = 0; t_idx < samples->betas; t_idx++) {
            fprintf(output_file, "beta\n");
            fprintf(output_file, "%lf\n", samples->beta_vals[t_idx]);
            fprintf(output_file, "w_k,L_HS_mean,L_HS_std\n");
            for (int k = 0; k < samples->k_max; k++) {
                fprintf(output_file, "%lf,%lf,%lf\n", samples->w_k[t_idx][k], samples->L_HS_mean[t_idx][k], samples->L_HS_std[t_idx][k]);
            }
        }
#endif // L_HH
    // } else {
    //     printf("Error in opening the %s file. \n", file_name);
    //     exit(1);
    // }

    fclose(output_file);
    return file_name;
}
