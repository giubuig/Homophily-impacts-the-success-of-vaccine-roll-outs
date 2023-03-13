using LinearAlgebra, Random, StatsBase

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
