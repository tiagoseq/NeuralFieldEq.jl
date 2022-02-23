# probNFE has to receive Input1D/Input2D as inputs 
# probNFE throws the ProbOutput1D/ProbOutput2D as outputs (which is of the inputs of solveSNFE)

@doc raw"""
`Input1D` is a structure to wrap the inputs to pass to `probNFE` function.
The arguments are in the following order:
# Input1D arguments
- `α  :: AbstractFloat`
- `v  :: AbstractFloat`
- `V0 :: fV0`: Can be a number or a function
- `L  :: Number`
- `N  :: Integer`
- `T  :: AbstractFloat`
- `n  :: Integer`
- `I  :: fI`: External input function
- `K  :: fK`: Connectivity function
- `S  :: fS`: Firing rate function

Declaration of `I`, `K` and `S` functions:
- The function I must have the arguments: `I(x,t)`
- The function K must have the arguments: `K(x)`
- The function S must have the arguments: `S(V)`
"""
struct Input1D{T<:Real,P<:Signed,O,fI,fK,fS,fV0} <: OneDim
          α    :: T
          v    :: T
          V0   :: fV0
          L    :: O
          N    :: P
          T    :: T
          n    :: P
    extInput   :: fI
    kernel     :: fK
    firingRate :: fS
end

@doc raw"""
`Input2D` is a structure to wrap the inputs to pass to `probNFE` function.
# Input2D arguments
- `α  :: AbstractFloat`
- `v  :: AbstractFloat`
- `V0 :: fV0`: Can be a number or a function
- `L  :: Number`
- `N  :: Integer`
- `T  :: AbstractFloat`
- `n  :: Integer`
- `I  :: fI`: External input function
- `K  :: fK`: Connectivity function
- `S  :: fS`: Firing rate function

Declaration of `I`, `K` and `S` functions:
- The function I must have the arguments: `I(x,y,t)`
- The function K must have the arguments: `K(x,y)`
- The function S must have the arguments: `S(V)`
"""
struct Input2D{T<:Real,P<:Signed,O,fI,fK,fS,fV0} <: TwoDim
    α          :: T
    v          :: T
    V0         :: fV0
    L          :: O
    N          :: P
    T          :: T
    n          :: P
    extInput   :: fI
    kernel     :: fK
    firingRate :: fS
end


@doc raw"""
In the 1D case, the output of `probNFE` function is `ProbOutput1D` structure.
This structure contains data to pass as input to function `solveNFE`.
"""
struct ProbOutput1D{T<:AbstractFloat,P<:Vector{<:Complex},A<:Domain,fI,fS} <: OneDim
    Plan    :: FFTW.rFFTWPlan{Float64,-1,false,1,UnitRange{Int64}}
    PlanInv :: AbstractFFTs.ScaledPlan{Complex{Float64},FFTW.rFFTWPlan{Complex{Float64},1,false,1,UnitRange{Int64}},Float64}
    krings  :: P
    V0      :: Vector{T}
    v0      :: P
    sv      :: P
    Ω       :: A
    α       :: T
    I       :: fI
    S       :: fS
end

@doc raw"""
In the 2D case, the output of `probNFE` function is `ProbOutput2D` structure.
This structure contains data to pass as input to function `solveNFE`.
"""
struct ProbOutput2D{T<:AbstractFloat,P<:Matrix{<:Complex},A<:Domain,fI,fS} <: TwoDim
    Plan    :: FFTW.rFFTWPlan{Float64,-1,false,2,UnitRange{Int64}}
    PlanInv :: AbstractFFTs.ScaledPlan{Complex{Float64},FFTW.rFFTWPlan{Complex{Float64},1,false,2,UnitRange{Int64}},Float64}
    krings  :: P
    V0      :: Matrix{T}
    v0      :: P
    sv      :: P
    Ω       :: A
    α       :: T
    I       :: fI
    S       :: fS
end