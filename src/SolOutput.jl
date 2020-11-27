struct SolOutput{T<:Vector{Float64},P<:Matrix{Float64}}
    V::Array{Float64,3}
    meanV::P
    x::T
    y::T
    t::T
    tsaved::T
end

# Take BuildSolution as a callable function with 4 methods
function (a::SolOutput)(ti::AbstractFloat,p::Int) #ISTO NÃO ESTÁ BEM. não é t-1 mas sim t[i-1]!
    t = findall(x->x==ti,a.tsaved)
    N = size(a.V)[2]
    return a.V[(t-1)*N+1:N*t,:,p]
end
(a::SolOutput)(ti::Int,p::Int) = a.V[(ti-1)*size(a.V)[2]+1:size(a.V)[2]*ti,:,p]

function (b::SolOutput)(ti::AbstractFloat) #ISTO NÃO ESTÁ BEM. não é t-1 mas sim t[i-1]!
    t = findall(x->x==ti,b.tsaved)
    N = size(b.meanV)[2]
    return b.meanV[(t-1)*N+1:N*t,:]
end
(b::SolOutput)(ti::Int64) = b.meanV[(ti-1)*size(b.meanV)[2]+1:size(b.meanV)[2]*ti,:]