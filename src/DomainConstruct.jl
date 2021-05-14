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

    # constructor
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