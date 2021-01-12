struct Domain{T<:Real,P<:Real} <: TwoDim
    # Space parameters
    L::T
    N::T
    dx::Float64
    dy::Float64
    x::Vector{Float64}
    y::Vector{Float64}
    # Time parameters
    tspan::P
    n::T
    dt::Float64
    t::Vector{Float64}

    # constructor
    function Domain{T,P}(L::T,N::T,tspan::P,n::T) where {T<:Real,P<:Real}
        # Space discretisation
        x_inf = -L/2 # Inferior lim in x
        x_sup = L/2  # Superior lim in x
        dx = (x_sup-x_inf)/N # x step
        dy = dx    # y step
    
        # Mesh grid
        x = collect(x_inf:dx:x_sup-dx)
        y = x

        # Time discretisation
        dt = tspan/n             # Time step
        t  = collect(0:dt:tspan) # Time discretized
        
        new(L,N,dx,dy,x,y,tspan,n,dt,t)
    end
end
