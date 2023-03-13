# Homophily-impacts-the-success-of-vaccine-roll-outs

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6108284.svg)](https://doi.org/10.5281/zenodo.6108284)

'MeanField.jl' allows to reproduce all the mean-field results. Those results used in the plots are provided in 'plots'.

'VaccDist_algorithm.jl' contains the algorithm to distribute the vaccine over the network. 'getCoverage.jl' contains the code to generate and save an arbitrary number of vaccine distributions for given values of coverage and mixing rate. The distributions are jointly stored in a jld file. Once a jld file with the vaccine distributions has been created, Monte Carlo simulations can be run via the 'MonteCarlo.jl'.

All the real data used are in 'data'. In particular, the code for measuring the time evolution of the vaccination homophily for each country is in 'HomophilyFromData.R'.
