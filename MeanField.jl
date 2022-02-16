using DifferentialEquations

function getR(α, V, R0, ϵ)
    kNN = 1-α*V
    kNV = α*V
    kVN = α*(1-V)
    kVV = 1-α*(1-V)

    a = 1
    b = -(R0*(kNN+(1-ϵ)*kVV) - 2)
    c = ((1-ϵ)*(kNN*kVV-kVN*kNV)*R0^2 - R0*(kNN+(1-ϵ)*kVV) + 1)

    return 1 + 0.5*(-b+sqrt(b*b-4*a*c))/a
end

function getR_expl(α, V, R0, ϵ)   #explicit expression for R (eq. 11)
    return (R0/2)*(2 - ϵ - α*(1-ϵ*(1-V)) + sqrt(α^2*(1-ϵ*(1-V))^2 + 2*α*ϵ*((1-ϵ)*(1-V)-V) + ϵ^2))
end

function sirRun!(du,u,p,t)
    R0, ϵ, kVN, kNV, kVV, kNN = p
    IV, IN, SV, SN, RV, RN = u

    du[1] = dIV = - IV + R0*(1-ϵ)*SV*(kVV*IV + kVN*IN)
    du[2] = dIN = - IN + R0*SN*(kNV*IV + kNN*IN)
    du[3] = dSV = - R0*(1-ϵ)*SV*(kVV*IV + kVN*IN)
    du[4] = dSN = - R0*SN*(kNV*IV + kNN*IN)
    du[5] = dRV  = IV
    du[6] = dRN  = IN
end

function get_prev(α, V, R0, ϵ)
    kNN = 1-α*V
    kNV = α*V
    kVN = α*(1-V)
    kVV = 1-α*(1-V)

    p = [R0, ϵ, kVN, kNV, kVV, kNN]
    p0 = 1e-3
    u0 = [p0, p0, 1-p0, 1-p0, 0., 0.]
    tmax = 1000.
    tspan = (0.0,tmax)
    steps = 10000
    prob = ODEProblem(sirRun!,u0,tspan,p)

    sol = solve(prob, saveat = tmax/steps)
    return sol
end

function getPeak(α, V, R0, ϵ)
    res = get_prev(α, V, R0, ϵ)

    IV = res[1,:]
    IN = res[2,:]
    I = IV*V + IN*(1-V)

    peakIV = maximum(IV)
    peakIN = maximum(IN)
    peakI = maximum(I)

    l = length(I)
    ΔIV = res[1,2:l] - res[1,1:(l-1)]
    ΔIN = res[2,2:l] - res[2,1:(l-1)]
    ΔI = ΔIV*V + ΔIN*(1-V)

    peakIV_inc = maximum(ΔIV)
    peakIN_inc = maximum(ΔIN)
    peakI_inc = maximum(ΔI)

    return(peakIV, peakIN, peakI, peakIV_inc, peakIN_inc, peakI_inc)
end
