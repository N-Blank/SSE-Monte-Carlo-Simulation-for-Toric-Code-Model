import numpy as np
import sys
from itertools import product

SAVE_DIR = "tmp/"
FILENAME = str(sys.argv[1])
N_LEGS = 12

# Read the input file to get the system parameters
with open(FILENAME, "r") as f:
    line = f.readline().strip().split(", ")

L = int(line[0])
N = 2 * L ** 2

S = float(line[1])
Sz = [m for m in np.arange(- S, S + 1, 1)]

SAVE_FILENAME = "vtx_S" + str(S) + ".txt"
print(SAVE_FILENAME)

def Sp_op(state):
    if state[1] + 1 == state[0]:
        return np.sqrt((S - state[1]) * (S + state[1] + 1))
    return 0.0

def Sm_op(state):
    if state[1] - 1 == state[0]:
        return np.sqrt((S + state[1]) * (S - state[1] + 1))
    return 0.0
    
def Sz_op(state):
    if state[0] == state[1]:
        return 2 * state[0]
    return 0.0

def H(state1, state2):
    # < state1 | H | state2 >
    # Hb = As + Bp
    # vtx type: 
    #   0 diagonal
    #   1 - 16 off_diagonal
  
    
    term = 0.0
    allowed = False
    vtx_type = 0
    
    spin_i = (state1[0], state2[0])
    spin_j = (state1[1], state2[1])
    spin_k = (state1[2], state2[2])
    spin_l = (state1[3], state2[3])
    spin_m = (state1[4], state2[4])
    spin_n = (state1[5], state2[5])
    
    term1 = (Sp_op(spin_i) * Sp_op(spin_j) * Sp_op(spin_k) * Sp_op(spin_l))
    if term1 != 0.0:
        vtx_type = 1
        allowed = True
        term = term1

    term2 = (Sp_op(spin_i) * Sp_op(spin_j) * Sp_op(spin_k) * Sm_op(spin_l))
    if term2 != 0.0:
        vtx_type = 2
        allowed = True
        term = term2

    term3 = (Sp_op(spin_i) * Sp_op(spin_j) * Sm_op(spin_k) * Sp_op(spin_l))
    if term3 != 0.0:
        vtx_type = 3
        allowed = True
        term = term3

    term4 = (Sp_op(spin_i) * Sp_op(spin_j) * Sm_op(spin_k) * Sm_op(spin_l))
    if term4 != 0.0:
        vtx_type = 4
        allowed = True
        term = term4

    term5 = (Sp_op(spin_i) * Sm_op(spin_j) * Sp_op(spin_k) * Sp_op(spin_l))
    if term5 != 0.0:
        vtx_type = 5
        allowed = True
        term = term5

    term6 = (Sp_op(spin_i) * Sm_op(spin_j) * Sp_op(spin_k) * Sm_op(spin_l))
    if term6 != 0.0:
        vtx_type = 6
        allowed = True
        term = term6

    term7 = (Sp_op(spin_i) * Sm_op(spin_j) * Sm_op(spin_k) * Sp_op(spin_l))
    if term7 != 0.0:
        vtx_type = 7
        allowed = True
        term = term7

    term8 = (Sp_op(spin_i) * Sm_op(spin_j) * Sm_op(spin_k) * Sm_op(spin_l))
    if term8 != 0.0:
        vtx_type = 8
        allowed = True
        term = term8

    term9 = (Sm_op(spin_i) * Sp_op(spin_j) * Sp_op(spin_k) * Sp_op(spin_l))
    if term9 != 0.0:
        vtx_type = 9
        allowed = True
        term = term9

    term10 = (Sm_op(spin_i) * Sp_op(spin_j) * Sp_op(spin_k) * Sm_op(spin_l))
    if term10 != 0.0:
        vtx_type = 10
        allowed = True
        term = term10

    term11 = (Sm_op(spin_i) * Sp_op(spin_j) * Sm_op(spin_k) * Sp_op(spin_l))
    if term11 != 0.0:
        vtx_type = 11
        allowed = True
        term = term11

    term12 = (Sm_op(spin_i) * Sp_op(spin_j) * Sm_op(spin_k) * Sm_op(spin_l))
    if term12 != 0.0:
        vtx_type = 12
        allowed = True
        term = term12

    term13 = (Sm_op(spin_i) * Sm_op(spin_j) * Sp_op(spin_k) * Sp_op(spin_l))
    if term13 != 0.0:
        vtx_type = 13
        allowed = True
        term = term13

    term14 = (Sm_op(spin_i) * Sm_op(spin_j) * Sp_op(spin_k) * Sm_op(spin_l))
    if term14 != 0.0:
        vtx_type = 14
        allowed = True
        term = term14

    term15 = (Sm_op(spin_i) * Sm_op(spin_j) * Sm_op(spin_k) * Sp_op(spin_l))
    if term15 != 0.0:
        vtx_type = 15
        allowed = True
        term = term15

    term16 = (Sm_op(spin_i) * Sm_op(spin_j) * Sm_op(spin_k) * Sm_op(spin_l))
    if term16 != 0.0:
        vtx_type = 16
        allowed = True
        term = term16

    term0 = 1 + (Sz_op(spin_m) * Sz_op(spin_n) * Sz_op(spin_k) * Sz_op(spin_l))
    if state1 == state2:
        vtx_type = 0
        allowed = True
        term = term0

    if (state1[4] != state2[4] or state1[5] != state2[5]):
        allowed = False
    
    return term, allowed, vtx_type

# Generate the vertex information
# Generate all of the possible vertices
all_vertex_state = list(product(Sz, repeat = N_LEGS))
allowed_vertex_state = list()
vtx = list()

# Select the allowed vertices
N_DIAGRAMS = 0
for vertex in all_vertex_state:
    term, allowed, vtx_type = H(vertex[6:], vertex[0:6])

    if allowed:
        allowed_vertex_state.append(list(vertex))
        
        vtx.append(dict())
        vtx[N_DIAGRAMS]["idx"] = N_DIAGRAMS
        vtx[N_DIAGRAMS]["type"] = vtx_type
        vtx[N_DIAGRAMS]["spin"] = list(vertex)
        vtx[N_DIAGRAMS]["H"] = np.around(term, decimals=10)
        vtx[N_DIAGRAMS]["new_vtx"] = [-1, -1]
        N_DIAGRAMS += 1

# Create transition table
for i in range(N_DIAGRAMS):
    new_state = vtx[i]["spin"].copy()
    for e in range(2):
        if e == 0:
            for x in range(4):
                new_state[x + 6] = -vtx[i]["spin"][x + 6]
        else:
            for x in range(4):
                new_state[x] = -vtx[i]["spin"][x]
                        
        if new_state in allowed_vertex_state:
            vtx[i]["new_vtx"][e] = allowed_vertex_state.index(new_state)

# Save the results to file
with open(SAVE_DIR + SAVE_FILENAME, "w") as f:
    f.write(str(N_DIAGRAMS) + "\n")
    f.write("\n")
    
    for i in range(N_DIAGRAMS):
        f.write(str(i) + "\n")
        f.write(str(vtx[i]["type"]) + "\n")
        f.write(str(vtx[i]["H"]) + "\n")
        
        for l in range(N_LEGS):
            f.write(str(int(2 * vtx[i]["spin"][l])) + " ")
        f.write("\n")

        for e in range(2):
            f.write(str(vtx[i]["new_vtx"][e]) + " ")
        f.write("\n")
    
        f.write("\n")
