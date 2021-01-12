### ### ### ### 1D problem:
# Deterministic 1D case
struct SolveOutDet1D{T<:Vector{Float64}} <: OneDim
    V::T
    x::T
    y::T
    t::T
    tsaved::T
end

# Stochastic 1D case
struct SolveOutSto1D{P<:Matrix{Float64},T<:Vector{Float64}} <: OneDim
    V::P
    meanV::T
    x::T
    y::T
    t::T
    tsaved::T
end

# Make struct Det 1D case as a callable function (2 methods)
function (a::SolveOutDet1D)(ti::AbstractFloat) #ISTO NÃO ESTÁ BEM. não é t-1 mas sim t[i-1]!
    t = findall(x->x==ti,a.tsaved)
    N = size(a.V)[2]
    return a.V[(t-1)*N+1:N*t,:]
end
(a::SolveOutDet1D)(ti::Int,p::Int) = a.V[(ti-1)*size(a.V)[2]+1:size(a.V)[2]*ti,:]

# Make struct Sto 1D case as a callable function (4 methods)
function (a::SolveOutSto1D)(ti::AbstractFloat,p::Int) #ISTO NÃO ESTÁ BEM. não é t-1 mas sim t[i-1]!
    t = findall(x->x==ti,a.tsaved)
    N = size(a.V)[2]
    return a.V[(t-1)*N+1:N*t,:,p]
end
(a::SolveOutSto1D)(ti::Int,p::Int) = a.V[(ti-1)*size(a.V)[2]+1:size(a.V)[2]*ti,:,p]

function (b::SolveOutSto1D)(ti::AbstractFloat) #ISTO NÃO ESTÁ BEM. não é t-1 mas sim t[i-1]!
    t = findall(x->x==ti,b.tsaved)
    N = size(b.meanV)[2]
    return b.meanV[(t-1)*N+1:N*t,:]
end
(b::SolveOutSto1D)(ti::Int64) = b.meanV[(ti-1)*size(b.meanV)[2]+1:size(b.meanV)[2]*ti,:]


### ### ### ### 2D problem:
# Deterministic 2D case
struct SolveOutDet2D{P<:Matrix{Float64},T<:Vector{Float64}} <: TwoDim
    V::P
    x::T
    y::T
    t::T
    tsaved::T
end

# Stochastic 2D case
struct SolveOutSto2D{M<:Array{Float64,3},P<:Matrix{Float64},T<:Vector{Float64}} <: TwoDim
    V::M
    meanV::P
    x::T
    y::T
    t::T
    tsaved::T
end

# Make struct Det 2D case as a callable function (2 methods)
function (a::SolveOutDet2D)(ti::AbstractFloat) #ISTO NÃO ESTÁ BEM. não é t-1 mas sim t[i-1]!
    t = findall(x->x==ti,a.tsaved)
    N = size(a.V)[2]
    return a.V[(t-1)*N+1:N*t,:]
end
(a::SolveOutDet2D)(ti::Int) = a.V[(ti-1)*size(a.V)[2]+1:size(a.V)[2]*ti,:]

# Make struct Sto 2D case as a callable function (4 methods)
function (a::SolveOutSto2D)(ti::AbstractFloat,p::Int) #ISTO NÃO ESTÁ BEM. não é t-1 mas sim t[i-1]!
    t = findall(x->x==ti,a.tsaved)
    N = size(a.V)[2]
    return a.V[(t-1)*N+1:N*t,:,p]
end
(a::SolveOutSto2D)(ti::Int,p::Int) = a.V[(ti-1)*size(a.V)[2]+1:size(a.V)[2]*ti,:,p]

function (b::SolveOutSto2D)(ti::AbstractFloat) #ISTO NÃO ESTÁ BEM. não é t-1 mas sim t[i-1]!
    t = findall(x->x==ti,b.tsaved)
    N = size(b.meanV)[2]
    return b.meanV[(t-1)*N+1:N*t,:]
end
(b::SolveOutSto2D)(ti::Int64) = b.meanV[(ti-1)*size(b.meanV)[2]+1:size(b.meanV)[2]*ti,:]