using DelimitedFiles, CSV, JLD
using StatsBase
using LightGraphs
using SparseArrays
using DataFrames

temporal_edges = DelimitedFiles.readdlm("data/t_edges_dbm74.dat")

weights_mtx, w_edgelist = AggregateTemporalEdgeList(temporal_edges)
N = size(weights_mtx,1)

th = 0.01
n_steps_max = 100000
V = 0.5

n_runs = 1000
n_alpha = 8
storeCov = zeros(n_runs, n_alpha, N)
αs = range(0.295, 0.995, length = n_alpha)

for j in 1:n_alpha
    n = 1
    while n <= n_runs
        α, Cov = getDist(w_edgelist, weights_mtx, V, th, αs[j], n_steps_max)
        if abs(α - αs[j]) < th
            storeCov[n,j,:] .= Cov
            n += 1
        end
    end
end

save("data/V05_dbm74.jld", "storeV05", storeCov)
