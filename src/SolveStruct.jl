#=The four possible outputs of `solveNFE` (1D-det/1D-sto and 2D-det/2D-sto) are defined
as structures, equipped with useful fields to be used in post-processing, such as
plotting, find solution values in specific points and time, know the maximum and minimum, etc.
The 1D solutions are vectors of size `N`, and 2D are matrices of size `N x N`.
Methods are defined to facilitate the access to solution values.=#

# 1D domain
@doc raw"""
The fields available in the case of the 1D deterministic solution:
- `V`: Vector in 1D, with the solution saved at the specified time instants
- `x`: Vector with the discretised space in direction $x$
- `tsaved`: Time instants where the solution was saved (for auxiliary purposes)
- `N`: Dimension of the solution at a specific time instant

## Methods available:
The domain was discretised as `x=x1,...,xN`

Solve the 1D NFE: `V = solveNFE(prob,[tj,tk,tl])`.
- `V(tj::Float)`: Returns the solution at `tj` (vector)
- `V(1::Int)`: Returns the solution at `tj` by referring to its index (V(tj)==V(1))
- `V(xi::Float,tk)`: Returns 1D solution at point `xi` at `tk`
"""
struct SolveOutDet1D{T<:Vector{<:AbstractFloat},C<:AbstractVector{<:AbstractFloat},P<:Signed} <: OneDim
    V      :: T
    x      :: T
    t      :: T
    tsaved :: C
    N      :: P
end

@doc raw"""
The fields available in the case of the 1D stochastic solution:
- `V`: Matrix `(N x np)`, where the `p-th` column corresponds to the `p-th` trajectory.
- `meanV`: Matrix with the mean solution across all `np` trajectories
- `x`: Vector with the discretised space in direction $x$
- `tsaved`: Time instants where the solution was saved (auxiliary purposes)
- `N`: Dimension of the solution at a specific time instant (auxiliary purposes)

## Methods available:
The domain was discretised as `x=x1,...,xN`

Solve the 1D SNFE: `Vsto = solveNFE(prob,[tj,tk,tl],ϵ,np)`.
- `Vsto(tj::Float,p::Int)`: Returns the `p-th` trajectory at `tj` (vector)
- `Vsto(tj::Float)`: Returns the mean solution at `tj` (vector)
- `Vsto(xi,tk,p)`: Returns `p-th` trajectory at point `xi` at `tk`
- `Vsto(xi,tk)`: Returns mean solution at point `xi` at `tk`
"""
struct SolveOutSto1D{T<:Vector{<:AbstractFloat},C<:AbstractVector{<:AbstractFloat},P<:Matrix{<:AbstractFloat},O<:Signed} <: OneDim
    V      :: P
    meanV  :: T
    x      :: T
    t      :: T
    tsaved :: C
    N      :: O
end

include("SolveStructMethods1D.jl")

# 2D domain
@doc raw"""
The fields available in the case of the 2D deterministic solution:
- `V`: Matrix with the solution saved at the specified time instants
- `x`: Vector with the discretised space in direction $x$
- `y`: Vector with the discretised space in direction $y$
- `tsaved`: Time instants where the solution was saved (auxiliary purposes)
- `N`: Matrix of size N x N (auxiliary purposes).

## Methods available:
The domain was discretised as `x=x1,...,xN` and `y=y1,...,yN`.

Solve the 2D NFE: `V = solveNFE(prob,[tj,tk,tl])`.
- `V(tj::Float)`: Returns the solution at `tj` (matrix)
- `V(1::Int)`: Returns the solution at `tj` by referring to its index (V(tj)==V(1))
- `V(xi::Float,yj::Float,tk)`: Returns solution at point `(xi,yj)` at `tk`
"""
struct SolveOutDet2D{P<:Matrix{<:AbstractFloat},T<:Vector{<:AbstractFloat},C<:AbstractVector{<:AbstractFloat},O<:Signed} <: TwoDim
    V      :: P
    x      :: T
    y      :: T
    t      :: T
    tsaved :: C
    N      :: O
end

@doc raw"""
The fields available in the case of the 2D stochastic solution:
- `V`: 3-dimensional array `(N x (N x tsaved) x np)`, where submatrix `N x N` corresponds
to the `p-th` trajectory at `tsaved` instant
- `meanV`: Matrix with the mean solution across all `np` trajectories
- `x`: Vector with the discretised space in direction $x$
- `y`: Vector with the discretised space in direction $y$
- `tsaved`: Time instants where the solution was saved (auxiliary purposes)
- `N`: Matrix of size N x N (auxiliary purposes).

## Methods available:
The domain was discretised as `x=x1,...,xN` and `y=y1,...,yN`.

Solve the 2D SNFE: `Vsto = solveNFE(prob,[tj,tk,tl],ϵ,np)`.
- `Vsto(tj::Float,p::Int)`: Returns the `p-th` trajectory at `tj` (matrix)
- `Vsto(tj::Float)`: Returns the mean solution at `tj` (matrix)
- `Vsto(xi,yj,tk,p)`: Returns `p-th` trajectory at point `(xi,yj)` at `tk`
- `Vsto(xi,yj,tk)`: Returns mean solution at point `(xi,yj)` at `tk`
"""
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