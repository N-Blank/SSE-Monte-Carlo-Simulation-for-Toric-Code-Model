#ifndef IO_H
#define IO_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "../sampling/sampling.h"
#include "../vtx/vtx_type.h"

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
void read_inputs(int *L, double *S, long *therm_cycles, long *mc_cycles, 
    int *n_bins, double **beta_vals, int *len_beta);

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
void read_vtx_info(char *file_name, vtx_element **vtx, int *n_diagrams);

/*
 * function: write_outputs
 *  writes the simulation outputs to file
 *   
 *  parameters:
 *      (sampled_quantities *) samples: sampled quantities
 *  during the simulation
 *  returns:
 *      (char *) file_name: name of save file
 */
char *write_outputs(sampled_quantities *samples, 
     int L, double S, long therm_cycles, long mc_cycles,
      double cpu_time_used, int n_threads);

#endif // IO_H
