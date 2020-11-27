struct ProbOutput{T<:AbstractFloat,P<:Int}
    P::FFTW.rFFTWPlan{Float64,-1,false,2,UnitRange{Int64}}
    Krings::Matrix{Complex{Float64}}
    rings::P
    I::Function
    S::Matrix{Complex{Float64}}
    firingRate::Function # The firing rate function
    Ω::Domain
    λ::Matrix{Float64}
    α::T
    ϵ::T
    np::P
    Vthresh::T
    V0::Matrix{Float64}
end