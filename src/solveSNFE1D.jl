# Deterministic method for 1D
function solveSNFE(prob::ProbOutput1D,saveat)
    P      = prob.Plan
    Pinv   = prob.PlanInv
    Krings = prob.Krings
    α      = prob.α
    rings  = prob.Ω.rings1D
    N      = prob.Ω.N
    t      = prob.Ω.t
    x      = prob.Ω.x
    dt     = prob.Ω.dt
    V      = prob.V0
    v      = prob.v0
    sv     = prob.sv # Firing rate history in Fourier Space (vector)
    hN     = N÷2+1   # Half of dim N (due to the output of rfft)

    # Pre-allocate vectors and matrices
    Vj  = Vector{Float64}(undef,N*length(saveat))
    svu = Vector{ComplexF64}(undef,hN) # Firing Rate in Fourier space (delay i)
    I   = similar(V)
    î   = similar(svu)
    a   = similar(svu) # Fourier coefficients of the integral operator A(X,t)

    j = 1 # Index for saveat
    @inbounds for i = 1:length(t) # Time loop
        ti = t[i]
        # Store solution V in saveat instants
        if ti in saveat
            Vj[(j-1)*N+1:N*j] .= V
            j += 1
        end

        I  = [prob.I(i,ti) for i in x]
        mul!(î,P,I)

        @. a = @views Krings[1:hN]*sv[1:hN]
        @inbounds for u = 2:rings
            @. a += @views Krings[1+hN*(u-1):hN*u]*sv[1+hN*(u-1):hN*u]
        end

        # Euler explicit scheme: v_i+1 = v_i + (Δt/α)*(î_i - v_i + a_i)
        @. v = v + (dt/α)*(î - v + a)

        V .= Pinv * v # Perform inverse Fourier transform to compute V in natural domain
        V .= V

        # Update s
        sv[hN+1:end] .= @view sv[1:hN*(rings-1)]
        mul!(svu,P,prob.S.(V)) # Apply Real Fourier Transform
        sv[1:hN] .= svu
    end #end time loop

    # Build our solution structure output
    sol = SolveOutDet1D(Vj,x,t,saveat,N)
    
    return sol

end

# Stochastic method for 1D
function solveSNFE(prob::ProbOutput1D,saveat,ϵ,np,ξ=0.1)
    P      = prob.Plan
    Pinv   = prob.PlanInv
    Krings = prob.Krings
    α      = prob.α
    rings  = prob.Ω.rings1D
    L      = prob.Ω.L
    N      = prob.Ω.N
    t      = prob.Ω.t
    x      = prob.Ω.x
    dt     = prob.Ω.dt
    hN     = N÷2+1   # Half of dim N (due to the output of rfft)
    
    # Pre-allocate matrices
    meanVj = Vector{Float64}(undef,N*length(saveat))
    Vj  = Matrix{Float64}(undef,N*length(saveat),np)
    V   = Vector{Float64}(undef,N)
    sv  = Vector{ComplexF64}(undef,hN*rings)
    svu = Vector{Complex{Float64}}(undef,hN) # Firing rate history in Fourier Space (vector)
    I   = similar(V)
    v   = similar(svu)
    î   = similar(svu)
    a   = similar(svu) # Fourier coefficients of the integral operator A(X,t)
    λ   = [exp(-(i^2*ξ^2)/(8*pi)) for i = 1:hN] # Correlation matrix
    
    # Pre-compute constant term of the stochastic part
    # fft gives the fourier coeff mult by N. Need to scale up noise, to have same magnitude
    # √dt due to the explicit Euler-Maruyama
    # Shift λ due to the also shifted output of fft
    c = N*sqrt(dt)*ϵ*fftshift(λ) # ϵ is the level of additive noise
    
    # Trajectories loop
    @inbounds for p = 1:np
        copy!(V,prob.V0) # Initial condition V(X,0) = V0 (Natural domain)
        copy!(v,prob.v0) # Initial condition v(x,0) = v0 (Fourier domain)
        copy!(sv,prob.sv)  # Initialise with the initial values of S
        
        j = 1 # Index for tj
        @inbounds for i = 1:length(t) # Time loop
            w = randn(hN) + im*randn(hN) # Stochastic term. w ~ N(0,1) + i*N(0,1) 
            ti = t[i]
            # Store solution V in saveat instants
            if ti in saveat
                Vj[(j-1)*N+1:N*j,p] .= V
                j += 1
            end
    
            I  = [prob.I(i,ti) for i in x]
            mul!(î,P,I)
    
            @. a = @views Krings[1:hN]*sv[1:hN]
            @inbounds for u = 2:rings
                @. a += @views Krings[1+hN*(u-1):hN*u]*sv[1+hN*(u-1):hN*u]
            end

            # Euler-Maruyama explicit scheme: v_i+1 = v_i + (Δt/α)*(î_i - v_i + a_i) + √(Δt)*ϵ*λ*w
            @. v = v + (dt/α)*(î - v + a) + c*w

            V .= Pinv * v # Perform inverse Fourier transform to compute V in natural domain

            # Update s
            sv[hN+1:end] .= @view sv[1:hN*(rings-1)]
            mul!(svu,P,prob.S.(V)) # Apply the Real Fourier Transform to V
            sv[1:hN,:].= svu       # Update s first delay ring
        end #end time loop

    end #end trajectories loop

    # Compute the mean solution across all trajectories
    @inbounds for i = 1:N*length(saveat)
        meanVj[i] = mean(@view Vj[i,:])
    end

    # Build our solution structure output
    sol = SolveOutSto1D(Vj,meanVj,x,t,saveat,N)
    
    return sol

end