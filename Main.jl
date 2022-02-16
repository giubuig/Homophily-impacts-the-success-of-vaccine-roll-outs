using LinearAlgebra, SparseArrays
using Random, Statistics, StatsBase
using DelimitedFiles, JLD

function montecarlo(edges_t, Vᵢ, β, μ, ϵ, p₀)
    N = length(Vᵢ)

    σ = ones(Int8,N)
    σ[sample(1:N, Int64(floor(p₀*N)), replace = false)] .= 2

    Itot = count(σ.==2)
    Rtot = count(σ.==3)
    Imax = copy(Itot)

    temp_σ = copy(σ)
    max_t = edges_t[end,1]
    L = size(edges_t, 1)

    t₀ = edges_t[1,1]
    tₑ = copy(t₀)
    indexₑ = 1
    βp = β*(1. - ϵ)

    while (Itot > 0)
        tₑ = copy(t₀)
        indexₑ = 1
        for t in t₀:max_t
            while (t == tₑ)
                e₁ = Int64(edges_t[indexₑ,2])
                e₂ = Int64(edges_t[indexₑ,3])
                if (σ[e₁] == 2 && temp_σ[e₂]*σ[e₂] == 1)
                    if Vᵢ[e₂] == 0.
                        if rand() < β
                            temp_σ[e₂] = 2
                        end
                    else
                        if rand() < βp
                            temp_σ[e₂] = 2
                        end
                    end
                elseif (temp_σ[e₁]*σ[e₁] == 1 && σ[e₂] == 2)
                    if Vᵢ[e₁] == 0.
                        if rand() < β
                            temp_σ[e₁] = 2
                        end
                    else
                        if rand() < βp
                            temp_σ[e₁] = 2
                        end
                    end
                end

                indexₑ += 1
                if indexₑ <= L
                    tₑ = edges_t[indexₑ,1]
                else
                    tₑ = -1
                end
            end

            @simd for i in 1:N
                if σ[i] == 2
                    if rand() < μ
                        temp_σ[i] = 3
                    end
                end
            end

            σ = copy(temp_σ)
            Itot = count(σ.==2)
            Rtot = count(σ.==3)
            if Itot > Imax
                Imax = copy(Itot)
            end

        end
    end

    return(Rtot, Imax)
end



#---  simulations

temporal_edges = DelimitedFiles.readdlm("data/t_edges_dbm74.dat")
secs_per_step = 300
temporal_edges[:,1] = temporal_edges[:,1]/secs_per_step .+ 1    #number of unique timestamps

weights_mtx, w_edgelist = AggregateTemporalEdgeList(temporal_edges)
N = size(weights_mtx,1)
strengths = [sum(weights_mtx[i,:]) for i in 1:N]
κ = mean((strengths./size(unique(temporal_edges[:,1]),1)).^2)/mean(strengths./size(unique(temporal_edges[:,1]),1))

R0 = 6.
μ = secs_per_step/(3600*24*7.5)
β = μ*R0/κ
ϵ = 0.6
# ϵ = 0.8
# ϵ = 1.0
p₀ = 10/N

V = 0.5
Cov = load("data/V05_dbm74.jld")["storeV05"]

n_runs = length(Cov[1:500,1,1])
n_alpha = length(Cov[1,:,1])
αs = range(0.3, 1.0, length = n_alpha)
n_runs_mc = 40

results = zeros(n_alpha, n_runs, n_runs_mc, 2)
res_inter = zeros(n_alpha, n_runs, 2)
resAttack_mc = zeros(Float64,length(αs),3,2)

@progress for indx in 1:length(αs)
    α_ = zeros(Float64,n_runs)
    @progress for r in 1:n_runs
        Vᵢ = Cov[r,indx,:]
        @progress for r_ in 1:n_runs_mc
            res = montecarlo(temporal_edges, Vᵢ, β, μ, ϵ, p₀)./N
            results[indx, r, r_, 1] = res[1]
            results[indx, r, r_, 2] = res[2]
        end
    end
end
for i in 1:n_alpha
    for j in 1:n_runs
        res_inter[i,j,1] = mean(results[i,j,:,1])
        res_inter[i,j,2] = mean(results[i,j,:,2])
    end
end
for i in 1:n_alpha
    for j in 1:2
        resAttack_mc[i,1,j] = median(res_inter[i,:,j])
        resAttack_mc[i,2,j] = quantile(res_inter[i,:,j], 0.25)
        resAttack_mc[i,3,j] = quantile(res_inter[i,:,j], 0.75)
    end
end
