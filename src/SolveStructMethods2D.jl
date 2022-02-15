# Struct Det 2D functor (4 methods) #

# Return V(x,y,ti)
function (a::SolveOutDet2D)(ti::AbstractFloat)
    @doc raw"""
    a(ti)

    This functor returns the 2D deterministic solution at instant `ti`
    """
    t = findlast(x->x==ti,a.tsaved)
    return a.V[(t-1)*a.N+1:a.N*t,:]
end
(a::SolveOutDet2D)(i::Int) = a.V[(i-1)*a.N+1:a.N*i,:]

# Return V(xj,yk,ti)
function (a::SolveOutDet2D)(xj::AbstractFloat,yk::AbstractFloat,ti::AbstractFloat)
    @doc raw"""
        a(xj,yk,ti)

    This functor returns the 2D deterministic solution at point (`xj`,`yk`) and instant `ti`
    """
    t  = findlast(i->i==ti,a.tsaved)
    x  = findlast(j->j==xj,a.x)
    y  = findlast(k->k==yk,a.y)
    Vi = @view a.V[(t-1)*a.N+1:a.N*t,:] # solution at ti
    
    return Vi[x,y]
end
function (a::SolveOutDet2D)(j::Int,k::Int,i::Int)
    Vi = @view a.V[(i-1)*a.N+1:a.N*i,:]
    return Vi[j,k]
end


# Struct Sto 2D functor (8 methods) #

# Return V(x,y,ti,p)
function (a::SolveOutSto2D)(ti::AbstractFloat,p::Int)
    @doc raw"""
        a(ti,p)

    This functor returns the path `p` at instant `ti` (2D domain)
    """
    t = findlast(i->i==ti,a.tsaved)
    return a.V[(t-1)*a.N+1:a.N*t,:,p]
end
(a::SolveOutSto2D)(i::Int,p::Int) = a.V[(i-1)*a.N+1:a.N*i,:,p]

# Return Vmean(x,y,ti)
function (b::SolveOutSto2D)(ti::AbstractFloat)
    @doc raw"""
        b(ti)

    This functor returns the 2D mean stochastic solution at instant `ti`
    """
    t = findlast(i->i==ti,b.tsaved)
    return b.meanV[(t-1)*b.N+1:b.N*t,:]
end
(b::SolveOutSto2D)(i::Int) = b.meanV[(i-1)*b.N+1:b.N*i,:]

# Return V(xj,yk,ti,p)
function (a::SolveOutSto2D)(xj::AbstractFloat,yk::AbstractFloat,ti::AbstractFloat,p::Int)
    @doc raw"""
        a(xj,yk,ti,p)

    This functor returns the path `p` at point (`xj`,`yk`) and instant `ti` (2D domain)
    """
    t  = findlast(i->i==ti,a.tsaved)
    x  = findlast(j->j==xj,a.x)
    y  = findlast(k->k==yk,a.y)
    Vi = @view a.V[(t-1)*a.N+1:a.N*t,:,p]
    
    return Vi[x,y]
end
function (a::SolveOutSto2D)(j::Int,k::Int,i::Int,p::Int)
    Vi = @view a.V[(i-1)*a.N+1:a.N*i,:,p]
    return Vi[j,k]
end

# Return Vmean(xj,yk,ti)
function (b::SolveOutSto2D)(xj::AbstractFloat,yk::AbstractFloat,ti::AbstractFloat)
    @doc raw"""
        b(xj,yk,ti)

    This functor returns the 2D mean stochastic solution at point (`xj`,`yk`) and instant `ti`
    """
    t  = findlast(i->i==ti,b.tsaved)
    x  = findlast(j->j==xj,b.x)
    y  = findlast(k->k==yk,b.y)
    Vi = @view b.meanV[(t-1)*b.N+1:b.N*t,:]
    
    return Vi[x,y]
end
function (b::SolveOutSto2D)(j::Int,k::Int,i::Int)
    Vi = @view b.meanV[(i-1)*b.N+1:b.N*i,:]
    return Vi[j,k]
end