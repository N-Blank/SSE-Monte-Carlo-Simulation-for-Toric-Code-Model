#include "sse.h"

/* 
 * function: diag_update
 *  diagonal update for the SSE MCS
 *  inserts or removes a diagonal operator from the operator string
 *  with some probablity
 * 
 *  paramters:
 *      (double) beta: temperature
 *      (toricode_system *) system: simulated system
 *      (sse_state *) state: SSE state
 */
void diag_update(double beta, toricode_system *system, sse_state *state, pcg32_random_t* rng) 
{
    for (int p = 0; p < state->M; p++) {
        if (state->op_string[p] == 0) {
            // no operator -> insert
            int b = pcg32_boundedrand_r(rng, system->Nb);
            if (pcg32_double_r(rng) <= (system->Nb * beta * prob(b, system, state)) / (state->M - state->n)) {
                state->op_string[p] = 17 * (b + 1);
                state->n++;
            }
        } else if (state->op_string[p] % 17 == 0) {
            // diagonal operator -> remove
            int b = (state->op_string[p] / 17) - 1;
            if (pcg32_double_r(rng) <= (state->M - state->n + 1) / (beta * system->Nb * 2)) {
                state->op_string[p] = 0;
                state->n--;
            }
        } else {
            // off-diagonal operator -> propagate

            int b = (state->op_string[p] / 17) - 1;

                system->spin[system->bond[b][0]] *= -1;
                system->spin[system->bond[b][1]] *= -1;
                system->spin[system->bond[b][2]] *= -1;
                system->spin[system->bond[b][3]] *= -1;

        }
    }
}


/* 
 * function: local off diagonal update
 *  local off diagonal update for the SSE MCS
 *  find the connected vertex pair and change them
 *  off-diagonal update
 * 
 *  parameters:
 *      (toricode_system *) system: simulated system
 *      (sse_state *) state: SSE state
 */
void loffdiagonal_update(toricode_system *system, sse_state *state, pcg32_random_t* rng) 
{
    // if no operators just randomize the spins
    if (state->n == 0) {
        for (int local = 0; local < state->n_locals; local++) {
            for (int i = 0; i < system->N; i++) {
                // update always but choose an update randomly
                system->spin[i] = system->Sz[pcg32_boundedrand_r(rng, system->n_proj)];
            }
        }
        return;
    }
    //bool flag = false;
    for (int p1 = 0; p1 < state->n; p1++) {
        for (int p2 = p1 + 1; p2 < state->n; p2++) {
            if (state->vtx_type[state->vtx[p1]].type == 0 && state->vtx_type[state->vtx[p2]].type != 0){
                continue;
            }
            if (state->vtx_type[state->vtx[p1]].type != 0 && state->vtx_type[state->vtx[p2]].type == 0){
                continue;
            }
            if (state->vtx_type[state->vtx[p1]].type == 0 && state->vtx_type[state->vtx[p2]].type == 0 && pcg32_double_r(rng) > 0.25){
                continue;
            }
            if (
                state->link[12 * p2] == 12 * p1 + 6
                &&state->link[12 * p2 + 1] == 12 * p1 + 7
                &&state->link[12 * p2 + 2] == 12 * p1 + 8
                &&state->link[12 * p2 + 3] == 12 * p1 + 9
                &&state->link[12 * p2 + 4] == 12 * p1 + 10
                &&state->link[12 * p2 + 5] == 12 * p1 + 11
            ){
                state->vtx[p1] = state->vtx_type[state->vtx[p1]].new_vtx_type[0];
                state->vtx[p2] = state->vtx_type[state->vtx[p2]].new_vtx_type[1];
                //printf("p1 %d \n", state->vtx_type[state->vtx[p1]].type);
                //printf("p2 %d \n", state->vtx_type[state->vtx[p2]].type);
                //flag = true;
                //break;
            }
            //if(flag)break;
        }
    }
    // translate the updated vertex list to the string operator
    for (int p = 0; p < state->n; p++) {
        state->red_op_string[p] = 17 * (state->red_op_string[p] / 17) + state->vtx_type[state->vtx[p]].type;
        state->op_string[state->trans_op_string[p]] = state->red_op_string[p];
    }

    // flip the update spins and flip randomly the non-updated ones
    for (int i = 0; i < system->N; i++) {
        if (state->first[i] != -1) {
            int p = state->first[i] / 12;
            int l = state->first[i] % 12;
            system->spin[i] = state->vtx_type[state->vtx[p]].spin[l];
        }
        else 
        {
            system->spin[i] = system->Sz[pcg32_boundedrand_r(rng, system->n_proj)];
        }
    }
}

/* 
 * function: ajust_cutoff
 *  dinamically adjusts the expansion cutoff during the thermalization
 *   part of the simulation
 * 
 *  parameters:
 *      (sse_state *) state: SSE state
 *      (long) t: mc time
 */
void ajust_cutoff(sse_state *state, long t) 
{
    u_int64_t M_new = state->n * 1.33;

    if (M_new > state->M) {
        u_int64_t *opstring_cpy = (u_int64_t *)malloc(state->M * sizeof(u_int64_t));
        memcpy(opstring_cpy, state->op_string, state->M * sizeof(u_int64_t));
        free(state->op_string); 
        state->op_string = (u_int64_t*) malloc(M_new * sizeof(u_int64_t));
        memset(state->op_string, 0, M_new * sizeof(u_int64_t));
        memcpy(state->op_string, opstring_cpy, state->M * sizeof(u_int64_t));
        free(opstring_cpy);

        state->M = M_new;
    }
}

/* 
 * function: init_toricode_system 
 *  initializes the toricode_system struct
 * 
 *  parameters:
 *      (int) L: number of unit cells
 *      (double) S: spin quantum number
 *      (toricode_system *) system: system to be initialized
 */
void init_toricode_system(int L, double S, toricode_system *system) 
{
    system->L = L;
    system->Nb = pow(L, 2);
    system->N = 2 * system->Nb;

    system->S = S;
    system->n_proj = 2.0 * S + 1.0;
    system->Sz = (int *) malloc(system->n_proj * sizeof(int));
    for (int i = 0; i < system->n_proj; i++) {
        system->Sz[i] = - 2 * (S - i);
    }
    system->spin = (int *) malloc(system->N * sizeof(int));
    system->bond = (int **) malloc(system->Nb * sizeof(int *));
    for (int i = 0; i < system->Nb; i++) { system->bond[i] = (int*) malloc(6 * sizeof(int)); }
    for (int i = 0; i < L; i++) {
        for (int j = 0; j < L; j++) {
            system->bond[j * L + i][4] = 2 * j * L + i;
            system->bond[j * L + i][5] = 2 * j * L + (i - 1 + L) % L;
            system->bond[j * L + i][0] = ((2 * j + 1) % (2 * L)) * L + i;
            system->bond[j * L + i][1] = ((2 * j - 1 + 2 * L ) % (2 * L)) * L + i;
            system->bond[j * L + i][2] = ((2 * (j + 1)) % (2 * L)) * L + i;
            system->bond[j * L + i][3] = ((2 * j + 1) % (2 * L)) * L + (i + 1) % L;
        }
    }
    for (int i = 0; i < system->N; i++) {
        system->spin[i] = system->Sz[0];
    }
}

/* 
 * function: init_sse_state 
 *  initializes the sse_state struct
 * 
 *  parameters:
 *      (uint64_t) seed: seed for the RNG
 *      (toricode_system *) system: simulated system
 *      (sse_state *) state: sse_state to be initialized
 */
void init_sse_state(uint64_t seed, toricode_system *system, sse_state *state) 
{
    state->n = 0;
    state->M = MAX_(12, system->N / 12);
    state->op_string = (u_int64_t *) malloc(state->M * sizeof(u_int64_t));
    memset(state->op_string, 0, state->M * sizeof(u_int64_t));

    state->first = (u_int64_t *) malloc(system->N * sizeof(u_int64_t));
    state->vtx = (u_int64_t *) malloc(state->n * sizeof(u_int64_t));
    state->link = (u_int64_t *) malloc(12 * state->n * sizeof(u_int64_t)); 
    state->red_op_string = (u_int64_t *) malloc(state->n * sizeof(u_int64_t));
    state->trans_op_string = (u_int64_t *) malloc(state->n * sizeof(u_int64_t));
}

/* 
 * function: reset_sse_state
 *  resets the SSE state
 *  deletes all of the operator string and spin states
 * 
 *  parameters:
 *      (toricode_system *) system: simulated system
 *      (sse_state *) state: SSE state
 */
void reset_sse_state(toricode_system *system, sse_state *state) 
{
    for (int i = 0; i < system->N; i++) {
        system->spin[i] = system->Sz[0];
    }

    state->n = 0;
    state->M = MAX_(12, system->N / 12);

    free(state->op_string);
    state->op_string = (u_int64_t *) malloc(state->M * sizeof(u_int64_t));
    memset(state->op_string, 0, state->M * sizeof(u_int64_t));
}

/* 
 * function: create_vtx_list
 *  tranlates the operator string into a doubly linked list of vertecies
 * 
 *  parameters:
 *      (toricode_system *) system: simulated system
 *      (sse_state *) state: SSE state
 */
void create_vtx_list(toricode_system *system, sse_state *state) 
{
    free(state->vtx);
    free(state->link);
    free(state->red_op_string);
    free(state->trans_op_string);
    state->vtx = (u_int64_t *) malloc(state->n * sizeof(u_int64_t));
    state->link = (u_int64_t *) malloc(12 * state->n * sizeof(u_int64_t)); 
    state->red_op_string = (u_int64_t *) malloc(state->n * sizeof(u_int64_t));
    state->trans_op_string = (u_int64_t *) malloc(state->n * sizeof(u_int64_t));
    
    u_int64_t last[system->N];
    memset(last, -1, system->N * sizeof(u_int64_t));
    memset(state->first, -1, system->N * sizeof(u_int64_t));
    memset(state->link, -1, 12 * state->n * sizeof(u_int64_t));
    memset(state->vtx, -1, state->n * sizeof(u_int64_t));
    memset(state->red_op_string, 0, state->n * sizeof(u_int64_t));
    memset(state->trans_op_string, 0, state->n * sizeof(u_int64_t));
    
    // reduce the operator string to the times where it has an operator
    // create a translation between the full and reduced one
    int p_red = 0;
    for (int p = 0; p < state->M; p++) {
        if (state->op_string[p] != 0) {
            state->red_op_string[p_red] = state->op_string[p];
            state->trans_op_string[p_red] = p;
            p_red++;
        }
    }

    // create the vertex list
    int l[12];
    for (int p = 0; p < state->n; p++) {
        // propagate the system
        int b = (state->red_op_string[p] / 17) - 1;
        int a = state->red_op_string[p] % 17;
        int vz = 12 * p;

        int i0 = system->bond[b][0];
        int i1 = system->bond[b][1];
        int i2 = system->bond[b][2];
        int i3 = system->bond[b][3];
        int i4 = system->bond[b][4];
        int i5 = system->bond[b][5];


        l[0] = system->spin[i0];
        l[1] = system->spin[i1];
        l[2] = system->spin[i2];
        l[3] = system->spin[i3];
        l[4] = system->spin[i4];
        l[5] = system->spin[i5];
        if (a != 0) {
            system->spin[i0] *= -1;
            system->spin[i1] *= -1;
            system->spin[i2] *= -1;
            system->spin[i3] *= -1;
        }
        l[6] = system->spin[i0];
        l[7] = system->spin[i1];
        l[8] = system->spin[i2];
        l[9] = system->spin[i3];
        l[10] = system->spin[i4];
        l[11] = system->spin[i5];

        // vertex at time p
        for (int i = 0; i < state->n_diagrams; i++) {
            if (l[0] == state->vtx_type[i].spin[0] && 
                l[1] == state->vtx_type[i].spin[1] && 
                l[2] == state->vtx_type[i].spin[2] && 
                l[3] == state->vtx_type[i].spin[3] && 
                l[4] == state->vtx_type[i].spin[4] && 
                l[5] == state->vtx_type[i].spin[5] && 
                l[6] == state->vtx_type[i].spin[6] && 
                l[7] == state->vtx_type[i].spin[7] && 
                l[8] == state->vtx_type[i].spin[8] && 
                l[9] == state->vtx_type[i].spin[9] && 
                l[10] == state->vtx_type[i].spin[10] && 
                l[11] == state->vtx_type[i].spin[11]){ state->vtx[p] = state->vtx_type[i].indx; }
        }
        // link the legs of the vertex
        int v0 = last[i0];
        int v1 = last[i1];
        int v2 = last[i2];
        int v3 = last[i3];
        int v4 = last[i4];
        int v5 = last[i5];
        
        if (v0 != -1) {
            state->link[v0] = vz;
            state->link[vz] = v0;
        } else {
            state->first[i0] = vz;
        }

        if (v1 != -1) {
            state->link[v1] = vz + 1;
            state->link[vz + 1] = v1;
        } else {
            state->first[i1] = vz + 1;
        }

        if (v2 != -1) {
            state->link[v2] = vz + 2;
            state->link[vz + 2] = v2;
        } else {
            state->first[i2] = vz + 2;
        }

        if (v3 != -1) {
            state->link[v3] = vz + 3;
            state->link[vz + 3] = v3;
        } else {
            state->first[i3] = vz + 3;
        }

        if (v4 != -1) {
            state->link[v4] = vz + 4;
            state->link[vz + 4] = v4;
        } else {
            state->first[i4] = vz + 4;
        }

        if (v5 != -1) {
            state->link[v5] = vz + 5;
            state->link[vz + 5] = v5;
        } else {
            state->first[i5] = vz + 5;
        }

        last[i0] = vz + 6;
        last[i1] = vz + 7;
        last[i2] = vz + 8;
        last[i3] = vz + 9;
        last[i4] = vz + 10;
        last[i5] = vz + 11;
    }
    
    for (int i = 0; i < system->N; i++) {
        if (state->first[i] != -1) {
            state->link[state->first[i]] = last[i];
            state->link[last[i]] = state->first[i];
        }
    }
}

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
double prob(int b, toricode_system *system, sse_state *state) 
{
    for (int i = 0; i < state->n_diagrams; i++) {
        if (state->vtx_type[i].type == 0 
        && state->vtx_type[i].spin[0] == system->spin[system->bond[b][0]] 
        && state->vtx_type[i].spin[1] == system->spin[system->bond[b][1]]
        && state->vtx_type[i].spin[2] == system->spin[system->bond[b][2]]
        && state->vtx_type[i].spin[3] == system->spin[system->bond[b][3]]
        && state->vtx_type[i].spin[4] == system->spin[system->bond[b][4]]
        && state->vtx_type[i].spin[5] == system->spin[system->bond[b][5]])
        {
            return state->vtx_type[i].H;
        }
    }
    return 0.0;
}

/* 
 * function: free_memory
 *  frees allocated memory in the toricode_system and sse_state structs
 * 
 *  parameters:
 *      (toricode_system *) system: toricode_system struct
 *      (sse_state *) state: sse_state struct
 */
void free_memory(toricode_system *system, sse_state *state) 
{
    for (int i = 0; i < system->Nb; i++) { free(system->bond[i]); }
    free(system->bond);
    free(system->spin);
    free(system->Sz);

    free(state->first);
    free(state->op_string);

    free(state->vtx);
    free(state->link);
    free(state->red_op_string);
    free(state->trans_op_string);

    free(state->vtx_type);
}
