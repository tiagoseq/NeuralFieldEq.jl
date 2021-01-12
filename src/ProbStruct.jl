# Structures for the probSNFE input
struct Input1D{T<:AbstractFloat,P<:Signed,fI,fK,fS} <: OneDim
    α::T
    v::T
    V0::T
    L::P
    N::P
    T::T
    n::P
    extInput::fI
    kernel::fK
    firingRate::fS
end

struct Input2D{T<:AbstractFloat,P<:Signed,fI,fK,fS} <: TwoDim
    α::T
    v::T
    V0::T
    L::P
    N::P
    T::T
    n::P
    extInput::fI
    kernel::fK
    firingRate::fS
end


# Structures for the probSNFE output
struct ProbOutput1D{P<:Signed} <: OneDim
    P::FFTW.rFFTWPlan{Float64,-1,false,1,UnitRange{Int64}}
    Krings::Matrix{Complex{Float64}}
    rings::P
    s::Vector{Complex{Float64}}
    V0::Vector{Float64}
    Ω::Domain
    in::Input1D
end

struct ProbOutput2D{P<:Signed} <: TwoDim
    P::FFTW.rFFTWPlan{Float64,-1,false,2,UnitRange{Int64}}
    Krings::Matrix{Complex{Float64}}
    rings::P
    s::Matrix{Complex{Float64}}
    V0::Matrix{Float64}
    Ω::Domain
    in::Input2D
end