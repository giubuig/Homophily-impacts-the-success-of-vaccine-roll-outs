# Homophily-impacts-the-success-of-vaccine-roll-outs

'MeanField.jl' allows to reproduce all the mean-field results.

The used temporal network is encoded as a time-ordered edgelist in the folder 'data'.

'VaccDist_algorithm.jl' contains the algorithm to distribute the vaccine over the network. 'getCoverage.jl' contains the code to generate and save an arbitrary number of vaccine distributions for given values of coverage and mixing rate. The distributions are jointly stored in a jld file. Once a jld file with the vaccine distributions has been created, Monte Carlo simulations can be run via the 'Main.jl'.
