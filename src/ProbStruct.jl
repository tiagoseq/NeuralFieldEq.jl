# Structures for the probSNFE input
struct Input1D{T<:Real,P<:Signed,fI,fK,fS} <: OneDim
    α::T
    v::T
    V0
    L::P
    N::P
    T::T
    n::P
    extInput::fI
    kernel::fK
    firingRate::fS
end

struct Input2D{T<:Real,P<:Signed,fI,fK,fS} <: TwoDim
    α::T
    v::T
    V0
    L::P
    N::P
    T::T
    n::P
    extInput::fI
    kernel::fK
    firingRate::fS
end


# Structures for the probSNFE output
struct ProbOutput1D{T<:AbstractFloat,P<:Vector{<:Complex},A<:Domain,fI,fS} <: OneDim
    Plan::FFTW.rFFTWPlan{Float64,-1,false,1,UnitRange{Int64}}
    PlanInv::AbstractFFTs.ScaledPlan{Complex{Float64},FFTW.rFFTWPlan{Complex{Float64},1,false,1,UnitRange{Int64}},Float64}
    Krings::P
    v0::P
    s::P
    Ω::A
    α::T
    I::fI
    S::fS
end

struct ProbOutput2D{T<:AbstractFloat,P<:Matrix{<:Complex},A<:Domain,fI,fS} <: TwoDim
    Plan::FFTW.rFFTWPlan{Float64,-1,false,2,UnitRange{Int64}}
    PlanInv::AbstractFFTs.ScaledPlan{Complex{Float64},FFTW.rFFTWPlan{Complex{Float64},1,false,2,UnitRange{Int64}},Float64}
    Krings::P
    v0::P
    s::P
    Ω::A
    α::T
    I::fI
    S::fS
end