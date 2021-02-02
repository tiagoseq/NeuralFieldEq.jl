# Deterministic method for 1D
function solveSNFE(prob::ProbOutput1D,saveat)
    P      = prob.Plan
    rings  = prob.Ω.rings
    Krings = prob.Krings
    α      = prob.α
    N      = prob.Ω.N
    t      = prob.Ω.t
    x      = prob.Ω.x
    dt     = prob.Ω.dt
    v      = prob.v0
    s      = prob.s # Firing rate history in Fourier Space (vector)
    hN     = N÷2+1 # Half of dim N (due to the output of rfft)

    # Pre-allocate vectors and matrices
    Vj = Vector{Float64}(undef,N*length(saveat))
    V  = Vector{Float64}(undef,N)
    su = Vector{ComplexF64}(undef,hN) # Firing Rate in Fourier space (delay i)
    I  = similar(V)
    dv = similar(su)
    î  = similar(su)
    a  = similar(su) # kernel by firing rate product in Fourier space

    j = 0 # Index for saveat
    @inbounds for i = 1:length(t) # Time loop
        ti = t[i]
        # Store solution V in saveat instants
        if t[i] in saveat
            j += 1
            Vj[(j-1)*N+1:N*j] .= V
        end

        I  = [prob.I(i,ti) for i in x]
        mul!(î,P,ifftshift(I))
        î .= fftshift(î)

        @. a = @views Krings[1:hN]*s[1:hN]
        @inbounds for u = 2:rings
            @. a += @views Krings[1+hN*(u-1):hN*u]*s[1+hN*(u-1):hN*u]
        end

        # Euler explicit scheme: V_i+1 = V_i + (Δt/α)*(I_i - V_i + A_i)
        @. v = v + (dt/α)*(î - v + a)

        ldiv!(V,P,ifftshift(v))
        V .= fftshift(V) # Perform inverse Fourier transform to compute V in natural domain

        # Update s
        s[hN+1:end] .= @view s[1:hN*(rings-1)]
        mul!(su,P,ifftshift(prob.S.(V))) # Apply Real Fourier Transform
        s[1:hN] .= fftshift(su)
    end #end time loop

    # Build our solution structure output
    sol = SolveOutDet1D(Vj,x,t,saveat,N)

    println("Solution saved at time instants: $(saveat)")
    
    return sol

end

# Stochastic method for 1D
function solveSNFE(prob::ProbOutput1D,saveat,ϵ,np,ξ=0.1)
    P      = prob.Plan
    rings  = prob.Ω.rings
    Krings = prob.Krings
    α      = prob.α
    N      = prob.Ω.N
    t      = prob.Ω.t
    x      = prob.Ω.x
    dt     = prob.Ω.dt
    v      = prob.v0
    s      = prob.s # Firing rate history in Fourier Space (vector)
    hN     = N÷2+1 # Half of dim N (due to the output of rfft)
    
    # Pre-allocate matrices
    meanVj = Vector{Float64}(undef,N*length(saveat))
    Vj = Matrix{Float64}(undef,N*length(saveat),np)
    V  = Matrix{Float64}(undef,N)
    su = Vector{Complex{Float64}}(undef,hN) # Firing Rate in Fourier space (delay i)
    dv = similar(su)
    î  = similar(su)
    a  = similar(su) # kernel by firing rate product in Fourier space
    I  = similar(V)
    λ  = [exp(-(i^2*ξ^2)/(8*pi)) for i in x] # Correlation matrix
    
    # Trajectories loop
    @inbounds for p = 1:np
        copy!(v,prob.v0) # Initial condition V(X,0) = V0
        copy!(s,prob.s)  # Initialise with the initial values of S
        
        j = 1 # Index for tj
        @inbounds for i = 1:length(t) # Time loop
            w = rand(Normal(),N) # Stochastic term
            ti = t[i]
            # Store solution V in saveat instants
            if t[i] in saveat
                j += 1
                Vj[(j-1)*N+1:N*j] .= V
            end
    
            I  = [prob.I(i,ti) for i in x]
            mul!(î,P,ifftshift(I))
            î .= fftshift(î)
    
            @. a = @views Krings[1:hN]*s[1:hN]
            @inbounds for u = 2:rings
                @. a += @views Krings[1+hN*(u-1):hN*u]*s[1+hN*(u-1):hN*j]
            end

            # Euler explicit scheme: V_i+1 = V_i + (Δt/α)*(I_i - V_i + A_i) + √(Δt)*ϵ*λ*w
            @. v = v + (dt/α)*(î - v + a) + sqrt(dt)*ϵ*λ*w

            ldiv!(V,P,ifftshift(v))
            V .= fftshift(V) # Perform inverse Fourier transform to compute V in natural domain

            # Update s
            s[hN+1:end] .= @view s[1:hN*(rings-1)]
            mul!(su,P,ifftshift(prob.S.(V))) # Apply Real Fourier Transform
            s[1:hN] .= fftshift(su)

        end #end time loop

    end #end trajectories loop

    # Compute the mean solution across all trajectories
    @inbounds for i = 1:N*length(saveat)
        meanVj[1+hN*(i-1):hN*i] = mean(Vj[1+hN*(i-1):hN*i,:])
    end

    # Build our solution structure output
    sol = SolveOutSto1D(Vj,meanVj,x,t,saveat,N)

    println("Solution saved at time instants: $(saveat)")
    
    
    return sol

end