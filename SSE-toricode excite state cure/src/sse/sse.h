#ifndef SSE_H
#define SSE_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include <stdint.h>

#include "../vtx/vtx_type.h"
#include "../rng/pcg_basic.h"

#define MAX_(a, b) ((a) > (b) ? (a) : (b))
#define MIN_(a, b) ((a) < (b) ? (a) : (b))


/* 
 * struct: toricode_system
 *  information about the simulated system
 *
 *  parameters
 *      (int) L: number of unit cells
 *      (int) N: number of particles (pow(L, d))
 *      (int) Nb: number of bonds 
 *      (double) S: spin quantum number
 *      (int) n_proj: number of spin projections in the z-direction
 *      (int *) Sz: spin values in the z-direction
 *      (int *) spin: spin state (length N)
 *      (int **) bond: lattice information (length Nb x 2) 
 */
typedef struct toricode_system
{
    int L;
    int N;
    int Nb;

    double S;
    int n_proj;
    int *Sz;

    int *spin;
    int **bond;
} toricode_system;

/* 
 * struct: sse_state
 *  information about the state of the SSE simulation
 * 
 *  parameters:
 *      (int *) op_string: operator string (length M) 
 *      (int) M: maximum expansion of the partiton function
 *      (int) n: number of operators in the operator string
 *      (int *) link: linked list (length 12*n)
 *      (int *) first: first imaginary time which the spins appear (length N)
 *      (int *) vtx: vertex type at each imaginery time (length n)
 *      (vtx_element *) vtx_type: information about each vertex type (length n_diagrams) 
 *      (int) n_diagrams: number of available vtx_elements
 *      (int *) red_op_string: reduced operator string (length n)
 *      (int *) trans_op_string: translation between the reduced and normal operator string (length n) 
 */
typedef struct sse_state 
{
    u_int64_t *op_string;
    u_int64_t M;
    u_int64_t n;

    u_int64_t n_locals;
    u_int64_t *link;
    u_int64_t *first;
    u_int64_t *vtx;
    vtx_element *vtx_type;
    int n_diagrams;

    u_int64_t* red_op_string;
    u_int64_t* trans_op_string;
} sse_state;

/* 
 * function: init_toricode_system
 *  initializes the toricode_system struct
 * 
 *  parameters:
 *      (int) L: number of unit cells
 *      (double) S: spin quantum number
 *      (toricode_system *) system: system to be initialized
 */
void init_toricode_system(int L, double S, toricode_system *system);

/* 
 * function: init_sse_state 
 *  initializes the sse_state struct
 * 
 *  parameters:
 *      (uint64_t) seed: seed for the RNG
 *      (heisenberg_system *) system: simulated system
 *      (sse_state *) state: sse_state to be initialized
 */
void init_sse_state(uint64_t seed, toricode_system *system, sse_state *state);

/* 
 * function: reset_sse_state
 *  resets the SSE state
 *  deletes all of the operator strings and spin states
 * 
 *  parameters:
 *      (toricode_system *) system: simulated system
 *      (sse_state *) state: SSE state
 */
void reset_sse_state(toricode_system *system, sse_state *state);

/* 
 * function: diag_update
 *  diagonal update for the SSE MCS
 *  inserts or removes a diagonal operator from the operator string
 *   with some probablity
 * 
 *  paramters:
 *      (double) beta: temperature
 *      (toricode_system *) system: simulated system
 *      (sse_state *) state: SSE state
 */
void diag_update(double beta, toricode_system *system, sse_state *state, pcg32_random_t* rng);

/* 
 * function: loffdiagonal_update
 *  local off diagonal update for the SSE MCS
 *  off-diagonal update
 * 
 *  parameters:
 *      (toricode_system *) system: simulated system
 *      (sse_state *) state: SSE state
 */
void loffdiagonal_update(toricode_system *system, sse_state *state, pcg32_random_t* rng); 

/* 
 * function: ajust_cutoff
 *  dynamically adjusts the expansion cutoff during the thermalization
 *   part of the simulation
 * 
 *  parameters:
 *      (sse_state *) state: SSE state
 *      (long) t: mc time
 */
void ajust_cutoff(sse_state *state, long t);

/* 
 * function: create_vtx_list
 *  tranlates the operator string into a doubly linked list of vertecies
 * 
 *  parameters:
 *      (toricode_system *) system: simulated system
 *      (sse_state *) state: SSE state
 */
void create_vtx_list(toricode_system *system, sse_state *state);

/* 
 * function: prob
 *  computes probability for diagonal updates
 *  
 *  parameters:
 *      (int) b: bond for the insertion/removal of operator
 *      (toricode_system *) system: system
 *      (sse_state* ) state: SSE state
 *  
 *  returns:
 *      (double) prob: probability for the update
 */
double prob(int b, toricode_system *system, sse_state *state);

/* 
 * function: free_memory
 *  frees allocated memory in the toricode_system and sse_state structs
 * 
 *  parameters:
 *      (toricode_system *) system: toricode_system struct
 *      (sse_state *) state: sse_state struct
 */
void free_memory(toricode_system *system, sse_state *state);

#endif // SSE_H
