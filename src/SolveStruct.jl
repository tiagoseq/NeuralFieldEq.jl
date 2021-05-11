### ### ### ### 1D problem:
# Deterministic 1D case
struct SolveOutDet1D{T<:Vector{<:AbstractFloat},P<:Signed} <: OneDim
    V::T
    x::T
    t::T
    tsaved::T
    N::P
end

# Stochastic 1D case
struct SolveOutSto1D{T<:Vector{<:AbstractFloat},P<:Matrix{<:AbstractFloat},O<:Signed} <: OneDim
    V::P
    meanV::T
    x::T
    t::T
    tsaved::T
    N::O
end

# Make struct Det 1D case as a callable function (2 methods)
function (a::SolveOutDet1D)(ti::AbstractFloat)
    t = findlast(x->x==ti,a.tsaved)
    return a.V[(t-1)*a.N+1:a.N*t]
end
(a::SolveOutDet1D)(ti::Signed) = a.V[(ti-1)*a.N+1:a.N*ti]

# Make struct Sto 1D case as a callable function (4 methods)
function (a::SolveOutSto1D)(ti::AbstractFloat,p::Int)
    t = findlast(x->x==ti,a.tsaved)
    return a.V[(t-1)*a.N+1:a.N*t,p]
end
(a::SolveOutSto1D)(ti::Int,p::Int) = a.V[(ti-1)*a.N+1:a.N*ti,p]

function (b::SolveOutSto1D)(ti::AbstractFloat)
    t = findlast(x->x==ti,b.tsaved)
    return b.meanV[(t-1)*b.N+1:b.N*t]
end
(b::SolveOutSto1D)(ti::Int64) = b.meanV[(ti-1)*b.N+1:b.N*ti]


### ### ### ### 2D problem:
# Deterministic 2D case
struct SolveOutDet2D{P<:Matrix{<:AbstractFloat},T<:Vector{<:AbstractFloat},O<:Signed} <: TwoDim
    V::P
    x::T
    y::T
    t::T
    tsaved::T
    N::O
end

# Stochastic 2D case
struct SolveOutSto2D{M<:Array{<:AbstractFloat,3},P<:Matrix{<:AbstractFloat},T<:Vector{<:AbstractFloat},O<:Signed} <: TwoDim
    V::M
    meanV::P
    x::T
    y::T
    t::T
    tsaved::T
    N::O
end

# Make struct Det 2D case as a callable function (2 methods)
function (a::SolveOutDet2D)(ti::AbstractFloat)
    t = findlast(x->x==ti,a.tsaved)
    return a.V[(t-1)*a.N+1:a.N*t,:]
end
(a::SolveOutDet2D)(ti::Signed) = a.V[(ti-1)*a.N+1:a.N*ti,:]

# Make struct Sto 2D case as a callable function (4 methods)
function (a::SolveOutSto2D)(ti::AbstractFloat,p::Int)
    t = findlast(x->x==ti,a.tsaved)
    return a.V[(t-1)*a.N+1:a.N*t,:,p]
end
(a::SolveOutSto2D)(ti::Int,p::Int) = a.V[(ti-1)*a.N+1:a.N*ti,:,p]

function (b::SolveOutSto2D)(ti::AbstractFloat)
    t = findlast(x->x==ti,b.tsaved)
    return b.meanV[(t-1)*b.N+1:b.N*t,:]
end
(b::SolveOutSto2D)(ti::Int64) = b.meanV[(ti-1)*b.N+1:b.N*ti,:]