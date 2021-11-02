# Define 4 structures as the output of solveSNFE 1D-det/1D-sto; 2D-det/2D-sto
# Make these structures callable by defining methods
# in order to return the solution at a specific time instant and path
# Det example: t = [1.5,3.0,4.5,6.0]
# V(2) or V(3.0) -> returns the solution at t[2]
# Sto example:
# V(2,55) or V(3.0,55) -> returns the solution at t[2] at trajectory 55
# V(2) or V(3.0) -> returns the mean solution at t[2]

### ### ### ### 1D problem:
# Deterministic 1D case
struct SolveOutDet1D{T<:Vector{<:AbstractFloat},C<:AbstractVector{<:AbstractFloat},P<:Signed} <: OneDim
    V      :: T
    x      :: T
    t      :: T
    tsaved :: C
    N      :: P
end

# Stochastic 1D case
struct SolveOutSto1D{T<:Vector{<:AbstractFloat},C<:AbstractVector{<:AbstractFloat},P<:Matrix{<:AbstractFloat},O<:Signed} <: OneDim
    V      :: P
    meanV  :: T
    x      :: T
    t      :: T
    tsaved :: C
    N      :: O
end

include("SolveStructMethods1D.jl")


### ### ### ### 2D problem:
# Deterministic 2D case
struct SolveOutDet2D{P<:Matrix{<:AbstractFloat},T<:Vector{<:AbstractFloat},C<:AbstractVector{<:AbstractFloat},O<:Signed} <: TwoDim
    V      :: P
    x      :: T
    y      :: T
    t      :: T
    tsaved :: C
    N      :: O
end

# Stochastic 2D case
struct SolveOutSto2D{M<:Array{<:AbstractFloat,3},P<:Matrix{<:AbstractFloat},T<:Vector{<:AbstractFloat},C<:AbstractVector{<:AbstractFloat},O<:Signed} <: TwoDim
    V      :: M
    meanV  :: P
    x      :: T
    y      :: T
    t      :: T
    tsaved :: C
    N      :: O
end

include("SolveStructMethods2D.jl")