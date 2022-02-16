using StatsBase
using LightGraphs
using SparseArrays
using DataFrames

function mixing_relation(w_edgelist, Cov)
    L = size(w_edgelist,1)

    n = 0
    for l in 1:L
        if Cov[Int64(w_edgelist[l,1])] == Cov[Int64(w_edgelist[l,2])]
            n += w_edgelist[l,3]
        end
    end
    return n
end

function AggregateTemporalEdgeList(temporal_edges)
    N = length(unique(temporal_edges[:,2:3]))
    L = size(temporal_edges, 1)
    W = sparse(zeros(Int64,N,N))

    for l in 1:L
        i = Int64(temporal_edges[l,2])
        j = Int64(temporal_edges[l,3])
        W[i,j] += 1
        W[j,i] += 1
    end

    WL = zeros(Int64, 0, 3)
    for i in 1:N
        for j in (i+1):N
            w = W[i,j]
            if w > 0
                WL = vcat(WL, [i, j, w]')
            end
        end
    end

    return(W, WL)
end

function getDist(w_edgelist, w_mtx, V, th, α_targ, n_steps_max)
    N = size(w_mtx,1)
    Cov = zeros(Int64, N)
    Cov[sample(1:N, Int64(floor(V*N)), replace = false)] .= 1
    vacs = findall(x -> x == 1, Cov)
    no_vacs = findall(x -> x == 0, Cov)

    h = mixing_relation(w_edgelist, Cov)

    neighbs = [findall(x -> x > 0, w_mtx[i,:]) for i in 1:N]
    strengths = [sum(w_mtx[i,:]) for i in 1:N]

    w = sum(w_edgelist[:,3])
    α = (1-h/w)/(2*V*(1-V))

    n = 0
    k = 1
    index_vac = 0
    h_k_vac = 0
    h_k_novac = 0
    j = 1
    index_novac = 0
    h_j_vac = 0
    h_j_novac = 0

    nVacs = length(vacs)
    sampleVacs = collect(1:nVacs)
    noVacs = length(no_vacs)
    sampleNoVacs = collect(1:noVacs)
    strVac = mean(strengths[findall(x->x == 1, Cov)])
    strNoVac = mean(strengths[findall(x->x == 0, Cov)])

    @inbounds while (n < n_steps_max && (α-α_targ) > th)
        k = j
        while k == j
            index_vac = sample(sampleVacs)
            index_novac = sample(sampleNoVacs)
            k = vacs[index_vac]
            j = no_vacs[index_novac]
        end

        if (strVac >= strNoVac &&  strengths[j] <= strengths[k]) || (strVac <= strNoVac &&  strengths[j] >= strengths[k])

            h_k_vac = 0
            h_k_novac = 0
            for i in neighbs[k]
                w_ki = w_mtx[k,i]
                h_k_vac += Cov[i]*w_ki
                h_k_novac += (1 - Cov[i])*w_ki
                if i == j
                    h_k_novac += -1*w_ki
                end
            end
            h_j_vac = 0
            h_j_novac = 0
            for i in neighbs[j]
                w_ji = w_mtx[j,i]
                h_j_vac += Cov[i]*w_ji
                h_j_novac += (1 - Cov[i])*w_ji
                if i == k
                    h_j_vac += -1*w_ji
                end
            end
            if (h_j_vac + h_k_novac) > (h_k_vac + h_j_novac)
                Cov[j] = 1
                Cov[k] = 0

                vacs[index_vac] = j
                no_vacs[index_novac] = k

                strVac = strVac - (strengths[k] - strengths[j])/nVacs
                strNoVac = strNoVac - (strengths[j] - strengths[k])/noVacs

                h = h - (h_k_vac + h_j_novac) + (h_j_vac + h_k_novac)
                α = (1-h/w)/(2*V*(1-V))
            end
        end
        n += 1
    end
    return(α, Cov)
end
