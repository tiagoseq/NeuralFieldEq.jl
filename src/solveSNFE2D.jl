# Deterministic method for 2D
function solveSNFE(prob::ProbOutput2D,saveat)
    P      = prob.Plan
    rings  = prob.Ω.rings
    Krings = prob.Krings
    t      = prob.Ω.t
    x      = prob.Ω.x
    y      = prob.Ω.y
    α      = prob.α
    N      = prob.Ω.N
    dt     = prob.Ω.dt
    v      = prob.v0
    s      = prob.s # Firing rate history in Fourier Space (matrix)
    hN     = N÷2+1  # Half of dim N (due to the output of rfft)
    
    # Pre-allocate matrices
    Vj = Matrix{Float64}(undef,N*length(saveat),N) # Stores the solution at 'saveat' instants
    V  = Matrix{Float64}(undef,N,N)
    su = Matrix{ComplexF64}(undef,hN,N) # Fourier coefficients of the Firing Rate at delay u
    I  = similar(V)
    a  = similar(su) # Fourier coefficients of the integral operator A(X,t)
    î  = similar(su) # Fourier coefficients of the external input
    
    j = 1 # Index for saveat
    @inbounds for i = 1:length(t) # Time loop
        ti = t[i]
        # Store solution V in saveat instants
        if ti in saveat
            Vj[(j-1)*N+1:N*j,:] .= V
            j += 1
        end

        I = [prob.I(i,j,ti) for j in y, i in x]
        mul!(î,P,ifftshift(I))
        î .= fftshift(î)

        
        @. a = @views Krings[1:hN,:]*s[1:hN,:]
        @inbounds for u = 2:rings
            @. a += @views Krings[1+hN*(u-1):hN*u,:]*s[1+hN*(u-1):hN*u,:]
        end
        
        # Euler explicit scheme: V_i+1 = V_i + (Δt/α)*(I_i - V_i + A_i)
        @. v = v + (dt/α)*(î - v + a)

        ldiv!(V,P,ifftshift(v))
        V .= fftshift(V) # Perform inverse Fourier transform to compute V in natural domain

        s[hN+1:end,:].= @view s[1:hN*(rings-1),:] # Update last delay ring
        mul!(su,P,ifftshift(prob.S.(V))) # Apply the Real Fourier Transform to V
        s[1:hN,:].= fftshift(su) # Update s first delay ring

    end #end time loop

    # Build our solution structure output
    sol = SolveOutDet2D(Vj,prob.Ω.x,prob.Ω.y,prob.Ω.t,saveat)

    println("Solution saved at time instants: $(saveat)")
    
    return sol

end

# Stochastic method for 2D
function solveSNFE(prob::ProbOutput2D,saveat,ϵ,np,ξ=0.1)
    P      = prob.Plan
    rings  = prob.Ω.rings
    Krings = prob.Krings
    t      = prob.Ω.t
    x      = prob.Ω.x
    y      = prob.Ω.y
    α      = prob.α
    N      = prob.Ω.N
    dt     = prob.Ω.dt
    v      = prob.v0
    s      = prob.s # Firing rate history in Fourier Space (matrix)
    hN     = N÷2+1  # Half of dim N (due to the output of rfft)
    
    # Pre-allocate matrices
    meanVj = Matrix{Float64}(undef,N*length(saveat),N)
    Vj = Array{Float64,3}(undef,N*length(saveat),N,np) # Stores the solution at 'saveat' instants
    V  = Matrix{Float64}(undef,N,N)
    su = Matrix{ComplexF64}(undef,hN,N) # Fourier coefficients of the Firing Rate at delay u
    I  = similar(V)
    a  = similar(su) # Fourier coefficients of the integral operator A(X,t)
    î  = similar(su) # Fourier coefficients of the external input
    λ  = [exp(-((i^2+j^2)*ξ^2)/(8*pi)) for j in y, i in x] # Correlation matrix
    
    # Trajectories loop
    @inbounds for p = 1:np
        copy!(v,prob.v0) # Initial condition V(X,0) = V0
        copy!(s,prob.s)  # Initialise with the initial values of S

        j = 1 # Index for tj
        @inbounds for i = 1:length(t) # Time loop
            w = rand(Normal(),hN,N) # Stochastic term
            ti = t[i]
            # Store solution V in t[j] instant
            if t[i] in saveat
                Vj[(j-1)*N+1:N*j,:,p] .= V
                j += 1
            end

            I = [prob.I(i,j,ti) for j in y, i in x]
            mul!(î,P,ifftshift(I))
            î .= fftshift(î)

            @. a = @views Krings[1:hN,:]*s[1:hN,:]
            @inbounds for u = 2:rings
                @. a += @views Krings[1+hN*(u-1):hN*u,:]*s[1+hN*(u-1):hN*u,:]
            end

            # Euler explicit scheme: V_i+1 = V_i + (Δt/α)*(I_i - V_i + A_i) + √(Δt)*ϵ*λ*w
            @. v = v + (dt/α)*(î - v + a) + sqrt(dt)*ϵ*λ*w

            ldiv!(V,P,ifftshift(v))
            V .= fftshift(V) # Perform inverse Fourier transform to compute V in natural domain

            s[hN+1:end,:] .= @view s[1:hN*(rings-1),:] # Update last delay ring
            mul!(su,P,ifftshift(prob.S.(V))) # Apply the Real Fourier Transform to V
            s[1:hN,:] .= fftshift(su) # Update last delay ring

        end #end time loop

    end #end trajectories loop

    # Compute the mean solution across all trajectories
    @inbounds for j = 1:N
        @inbounds for i = 1:N*length(saveat)
            meanVj[i,j] = mean(Vj[i,j,:])
        end
    end

    # Build our solution structure output
    sol = SolveOutSto2D(Vj,meanVj,prob.Ω.x,prob.Ω.y,prob.Ω.t,saveat)

    println("Solution saved at time instants: $(saveat)")
    
    
    return sol

end