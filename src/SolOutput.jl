struct SolOutput
    V::Array{Float64,3}
    meanV::Array{Float64,2}
    x::Array{Float64,1}
    y::Array{Float64,1}
    t::Array{Float64,1}
    tsaved::Array{Float64,1}
end

# Take BuildSolution as a callable function with 4 methods
(a::SolOutput)(ti::Float64,p::Int64) = a.V[(findall(x->x==ti,a.t)-1)*size(a.V)[2]+1:size(a.V)[2]*findall(x->x==ti,a.t),:,p]
(a::SolOutput)(ti::Int64,p::Int64)   = a.V[(ti-1)*size(a.V)[2]+1:size(a.V)[2]*ti,:,p]
(b::SolOutput)(ti::Float64) = b.meanV[(findall(x->x==ti,b.t)-1)*size(b.meanV)[2]+1:size(b.meanV)[2]*findall(x->x==ti,b.t),:]
(b::SolOutput)(ti::Int64)   = b.meanV[(ti-1)*size(b.meanV)[2]+1:size(b.meanV)[2]*ti,:]