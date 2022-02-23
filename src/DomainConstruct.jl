@doc raw"""
This structure discretises the domain, in space and time,
and computes the discrete delay ring structure.
The fields of this structure are:
- `L` : Spatial domain length (`[-L/2 x L/2-dx]`)
- `N` : Number of spatial nodes, used for discretise space
- `dx`: Spatial step on `x` direction
- `dy`: Spatial step on `y` direction
- `x` : Discrete direction `x`
- `y` : Discrete direction `y`
- `T` : Time span of the simulation
- `n` : Number of time steps, used for discretise time
- `v` : Axonal velocity
- `Δr`: Discrete ring width
- `rings1D` : Number of delay rings in 1D domains
- `rings2D` : Number of delay rings in 2D domains

Containing all the data inherent to the discretised domain.
"""
struct Domain{T<:Real,P<:AbstractFloat,O<:Signed} <: TwoDim
    L       :: T
    N       :: O
    dx      :: P
    dy      :: P
    x       :: Vector
    y       :: Vector
    tspan   :: P
    n       :: O
    dt      :: P
    t       :: Vector
    v       :: P
    Δr      :: P
    rings1D :: O
    rings2D :: O

    # domain constructor
    function Domain{T,P,O}(L::T,N::O,tspan::P,n::O,v::P) where {T<:Real,P<:AbstractFloat,O<:Signed}
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

        # Delay computations
        τ_max_1D = L/(2*v)              # Maximum delay 1D domain
        τ_max_2D = L/(sqrt(2)*v)        # Maximum delay 2D domain
        Δr    = max(1.0,(v*dt)/dx)      # Compute rings width
        rings1D = 1+floor(O,τ_max_1D/dt) # Number of delay rings (u_max)
        rings2D = 1+floor(O,τ_max_2D/dt) # Number of delay rings (u_max)
        
        new(L,N,dx,dy,x,y,tspan,n,dt,t,v,Δr,rings1D,rings2D)
    end
end