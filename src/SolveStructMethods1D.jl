# Struct Det 1D functor (4 methods) #

# Return V(x,ti)
function (a::SolveOutDet1D)(ti::AbstractFloat)
    t = findlast(i->i==ti,a.tsaved)
    return a.V[(t-1)*a.N+1:a.N*t]
end
(a::SolveOutDet1D)(i::Int) = a.V[(i-1)*a.N+1:a.N*i]

# Return V(xj,ti)
function (a::SolveOutDet1D)(xj::AbstractFloat,ti::AbstractFloat)
    t  = findlast(i->i==ti,a.tsaved)
    x  = findlast(j->j==xj,a.x)
    Vi = @view a.V[(t-1)*a.N+1:a.N*t]
    
    return Vi[x]
end
function (a::SolveOutDet1D)(j::Int,i::Int)
    Vi = @view a.V[(i-1)*a.N+1:a.N*i] # solution at ti
    return Vi[j]
end


# Struct Sto 1D functor (8 methods) #

# Return V(x,ti,p)
function (a::SolveOutSto1D)(ti::AbstractFloat,p::Int)
    t = findlast(i->i==ti,a.tsaved)
    return a.V[(t-1)*a.N+1:a.N*t,p]
end
(a::SolveOutSto1D)(i::Int,p::Int) = a.V[(i-1)*a.N+1:a.N*i,p]

# Return Vmean(x,ti)
function (b::SolveOutSto1D)(ti::AbstractFloat)
    t = findlast(i->i==ti,b.tsaved)
    return b.meanV[(t-1)*b.N+1:b.N*t]
end
(b::SolveOutSto1D)(i::Int) = b.meanV[(i-1)*b.N+1:b.N*i]

# Return V(xj,ti,p)
function (a::SolveOutSto1D)(xj::AbstractFloat,ti::AbstractFloat,p::Int)
    t  = findlast(i->i==ti,a.tsaved)
    x  = findlast(j->j==xj,a.x)
    Vi = @view a.V[(t-1)*a.N+1:a.N*t,p] # path at ti
    
    return Vi[x] # solution at x and ti; i.e. V(x,t)
end
function (a::SolveOutSto1D)(j::Int,i::Int,p::Int)
    Vi = @view a.V[(i-1)*a.N+1:a.N*i,p]
    return Vi[j]
end

# Return Vmean(xj,ti)
function (b::SolveOutSto1D)(xj::AbstractFloat,ti::AbstractFloat)
    t  = findlast(i->i==ti,b.tsaved)
    x  = findlast(j->j==xj,b.x)
    Vi = @view b.meanV[(t-1)*b.N+1:b.N*t]
    
    return Vi[x]
end