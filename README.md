# SSE-Monte-Carlo-Simulation-for-Toric-Code-Model

The first folder `SEE-toricode` is exact implementation of Stochastic Series Expansion (SSE) Monte Carlo method for the toricode model. This version uses the Belt Loops [1] method for the loop update.
The Hamiltonian of the simulated system is given by

$$ H = - \sum_{\langle s \rangle} A_s - \sum_{\langle p \rangle} B_p $$

where $A_s = \prod_{j\in star(s)} \sigma_j^x$,  $B_p = \prod_{j\in boundary(p)} \sigma_j^z$, which is obviously the traditional toricode model.

The SSE method is the power series expansion of the partition function 

$$ Z = \text{Tr}\{e^{-\beta H}\} = \sum_{\alpha} \sum_{n=0}^{\infty} \frac{(-\beta)^n}{n!} \langle \alpha |H^n| \alpha \rangle $$

where $\beta \equiv 1 / T$.

The implementation uses binnig for the estimation of the standard deviations of the sampled quantities. It is possible to run each bin in parallel mode, using openMP.

## Usage

To use and run the implementation, it is required that you have a C compiler (gcc-12, preferably), and a version of Python3 installed. To install GCC and Python3 you run the following commands if you have a Debian based OS (Ubuntu, Pop_OS, ...)
```bash
$ sudo apt install gcc python3
```
or if you have a MacOS based device with [Homebrew](https://brew.sh)
```bash
$ brew install gcc python3
```
You will also need an uptodate version of NumPy. You can do so by running the following command
```bash
$ pip3 install numpy
```

To run a simulation you will need to use the `run.sh` script and three input files. The first input file `read.in` contains the information about the system and MC parameters for the simulation.
```
8, 0.5
10000, 1000000, 10
```
The first line is `L, S` and the second is `therm_cycles, mc_cycles, n_bins`. `L` is the number of unit cells, `S` is the quantum spin number.

The second input file is `beta.in` which has the following structure
```
8
0.5
1.0
2.0
4.0
```
The first line is the number of beta values in the file and the next lines are the beta values for the MC simulation.


After setting the simulation parameters, you will need to run it via the `run.sh` script. You can type 
```bash
$ ./run.sh -h
```
to get more information about the arguments. The most important ones are `-n <n_threads>` which specifies the number of threads used by openMP.
```bash
$ ./run.sh -n 5 
```
To see the plot of the exact result, run the python file `plot.py`.
```bash
$ python3 plot.py /L8_toricode_S0.5.csv
```
The second folder, titled "SEE-toricode excited state," contains the exact implementation of the SSE method for the excited toricode model. The model is described by the Hamiltonian:

$$ H = - \sum_{\langle s \rangle} A_s - \sum_{\langle p \rangle} B_p + 2A_v$$

Here, $A_v$ represents the $A_s$ operator at a specific vertex $v$. By studying this model, we can explore the extent to which the sign problem affects the results of the excited toricode model SSE simulation.

In the third folder, named "SEE-toricode excited state_cure," we demonstrate a method to address the sign problem in the Toric-code excited state by applying a basis change.
To observe the outcomes illustrated in the figures, you can employ the same approach as mentioned above to run the files.

[1] - "Stochastic Series Expansion Methods", Sandvik. Anders W., 2019,  http://arxiv.org/abs/1909.10591

[2] - "Stochastic Series Expansion Quantum Monte Carlo", Melko. Roger G., 2013, Springer Berlin Heidelberg, https://link.springer.com/10.1007/978-3-642-35106-8_7

[3] - "Fault-tolerant quantum computation by anyons", Kitaev,A. Yu., 2003, Annals of Physics, http://arxiv.org/abs/quant-ph/9707021

[4] - "Directed loop updates for quantum lattice models", Olav F. Syljuåsen, 2003, Phys. Rev. E 67, 046701, https://doi.org/10.1103/PhysRevE.67.046701

[5] - "Sign Problem in Quantum Monte Carlo Simulation", Gaopei Pan, Zi Yang Meng, 2022, https://arxiv.org/pdf/2204.08777.pdf
