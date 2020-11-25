struct ProbOutput
    P::FFTW.rFFTWPlan{Float64,-1,false,2,UnitRange{Int64}}
    Krings::Array{Complex{Float64},2}
    rings::Int64
    I::Function
    S::Array{Complex{Float64},2}
    firingRate::Function # The firing rate function
    Ω::Domain
    λ::Array{Float64,2}
    α::Float64
    ϵ::Float64
    np::Int64
    Vthresh::Float64
    V0::Array{Float64,2}
end