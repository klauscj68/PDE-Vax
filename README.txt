This repository contains the codes used in our PDE DSA paper applied to 
the COVID epidemic in Ohio from Nov 2020-Jan 2021. 

### Installing the environment
The Julia environment needed to run the code is setup using the Project and
Manifest files. One can activate the environment by switching to the REPL
package manager and using the activate and instantiate commands. See for 
example
https://pkgdocs.julialang.org/v1/environments/#Using-someone-else's-project.

Jupyter notebooks should automatically load installed pacakges without 
needing to activate the environment, although additional pacakges may need
to be instantiated.

### Running a single PDE DSA simulation
The simplest way to become familiar with the code is to run a simulation
from the jupyter notebook solve_pdedsa.ipynb. The notebook contains 
instructions on how to run the simulation. A sample output file 
solve_pdedsa_output.html is also provided.

Simulations may also be run in REPL by inputing the jupyter commands into
the julia command line. Alternative equation coefficients may be specified
by modifying the appropriate function in parameters.jl. 

### Running ABC
The julia files parameters.jl and abcmc.jl specify most of the options,
including fixed parameter values, priors, number of samples, and how the
error is computed. Fixed parameter values are taken from data() in
parameters.jl. Priors are specified in abcdata() and abcsmp!() in abcmc.jl.
The error is computed in ℓerr() in abcmc.jl. 

To generate the samples, one manually calls abcpregen() with the number of 
samples in a single batch and the number of batches. This is how 
simulations may be run embarassingly parallel. The code generates # batch 
many csv files for both sample parameter values and incidence trajectories.

Then to compute the ABC sample trajectories, one calls abcrun() with 
arguments specifying the specific sample and trajectory csv that the 
finished simulation should write to, which records the observed error and 
the incidence trajectory.

Once all # batch many runs have been run, the output is analyzed by the 
jupyter notebook process_abc.ipynb. A sample output file has been provided
in process_abc_output.html.

### Extracting parameter values from Ohio data
The notebook process_ODH.ipynb analyzes the Ohio epidemic data and extracts
parameteric and nonparametric values which were used in the ABC simulations.
Sample output has been provided in process_ODH_output.html.    
