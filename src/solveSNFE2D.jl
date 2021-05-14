# Deterministic method for 2D
function solveSNFE(prob::ProbOutput2D,saveat)
    P      = prob.Plan
    Pinv   = prob.PlanInv
    krings = prob.krings
    rings  = prob.Ω.rings2D
    t      = prob.Ω.t
    x      = prob.Ω.x
    y      = prob.Ω.y
    N      = prob.Ω.N
    dt     = prob.Ω.dt
    α      = prob.α
    V      = prob.V0 # Initial condition in natural space (needed if saveat contains t[1])
    v      = prob.v0 # Initial condition in Fourier space
    sv     = prob.sv # Firing Rate in Fourier space (all delays)
    hN     = N÷2+1   # Half of dim N (due to the output of rfft)
    
    # Pre-allocate matrices
    Vj  = Matrix{Float64}(undef,N*length(saveat),N) # Stores the solution at 'saveat' instants
    svu = Matrix{ComplexF64}(undef,hN,N) # Fourier coefficients of the Firing Rate at delay u
    I   = similar(V)
    î   = similar(svu) # Fourier coefficients of the external input
    a   = similar(svu) # Fourier coefficients of the integral operator A(X,t)
    j   = 1 # Set index for solution saved at saveat instants

    @inbounds for i = 1:length(t) # Time loop
        ti = t[i]
        # Store solution V in saveat instants
        if ti in saveat
            Vj[(j-1)*N+1:N*j,:] .= V
            j += 1
        end

        I.= [prob.I(i,j,ti) for j in y, i in x] # Discretise I
        mul!(î,P,I)                              # rfft of I

        # Compute the integral operator's fourier coef at t_i
        @. a = @views krings[1:hN,:]*sv[1:hN,:] # Init a
        @inbounds for u = 2:rings # Rings loop
            @. a += @views krings[1+hN*(u-1):hN*u,:]*sv[1+hN*(u-1):hN*u,:]
        end
        
        # Euler explicit scheme: v_i+1 = v_i + (Δt/α)*(î_i - v_i + a_i)
        @. v = v + (dt/α)*(î - v + a)

        V .= Pinv * v # Perform inverse Fourier transform to compute V in natural domain

        # Update sv
        # Every delay block moves one block down. The last one is deleted
        sv[hN+1:end,:].= @view sv[1:hN*(rings-1),:]
        mul!(svu,P,prob.S.(V)) # Apply the Real Fourier Transform to V
        sv[1:hN,:].= svu # Update sv first delay ring
    end #end time loop

    # Build our solution structure output
    sol = SolveOutDet2D(Vj,prob.Ω.x,prob.Ω.y,prob.Ω.t,saveat,N)
    
    return sol

end

# Stochastic method for 2D
function solveSNFE(prob::ProbOutput2D,saveat,ϵ,np,ξ=0.1)
    P      = prob.Plan
    Pinv   = prob.PlanInv
    krings = prob.krings
    rings  = prob.Ω.rings2D
    t      = prob.Ω.t
    x      = prob.Ω.x
    y      = prob.Ω.y
    N      = prob.Ω.N
    dt     = prob.Ω.dt
    α      = prob.α
    hN     = N÷2+1  # Half of dim N (due to the output of rfft)
    
    # Pre-allocate matrices
    meanVj = Matrix{Float64}(undef,N*length(saveat),N) # Store the mean solution at saveat
    Vj  = Array{Float64,3}(undef,N*length(saveat),N,np) # Store the solution at saveat at path p
    V   = Matrix{Float64}(undef,N,N)
    v   = Matrix{ComplexF64}(undef,hN,N)
    sv  = Matrix{ComplexF64}(undef,hN*rings,N) # Firing rate history in Fourier Space (matrix)
    svu = similar(v) # Fourier coefficients of the Firing Rate at delay u
    a   = similar(v) # Fourier coefficients of the integral operator A(X,t)
    î   = similar(v) # Fourier coefficients of the external input
    I   = similar(V)
    λ   = [exp(-((i^2+j^2)*ξ^2)/(8*pi)) for j = 1:hN, i = 1:N] # Correlation matrix

    # Pre-compute constant term of the stochastic part
    # fft gives the fourier coeff mult by N. Need to scale up noise, to have same magnitude
    # √dt due to the explicit Euler-Maruyama
    # Shift λ due to the also shifted output of fft
    c = (N^2)*sqrt(dt)*ϵ*fftshift(λ) # ϵ is the level of additive noise

    # Trajectories loop
    @inbounds for p = 1:np
        copy!(V,prob.V0) # Initial condition V(X,0) = V0 (Natural domain)
        copy!(v,prob.v0) # Initial condition v(X,0) = v0 (Fourier domain)
        copy!(sv,prob.sv)  # Initialise with the initial values of S

        j = 1 # Set index for tj
        @inbounds for i = 1:length(t) # Time loop
            # w ~ N(0,1) + i*N(0,1)  stochastic term
            w = randn(Float64,hN,N) + im*randn(Float64,hN,N)
            ti = t[i]
            # Store solution V in t[j] instant
            if ti in saveat
                Vj[(j-1)*N+1:N*j,:,p] .= V
                j += 1
            end

            I.= [prob.I(i,j,ti) for j in y, i in x] # Discretise I
            mul!(î,P,I)                             # rfft of I

            # Compute the integral operator's fourier coef at t_i
            @. a = @views krings[1:hN,:]*sv[1:hN,:] # Init a
            @inbounds for u = 2:rings # Ring loop
                @. a += @views krings[1+hN*(u-1):hN*u,:]*sv[1+hN*(u-1):hN*u,:]
            end

            # Euler explicit scheme: v_i+1 = v_i + (Δt/α)*(î_i - v_i + a_i) + √(Δt)*ϵ*λ*w
            @. v = v + (dt/α)*(î - v + a) + c*w
            
            # Apply inverse Fourier transform to compute V in natural domain
            V .= Pinv * v

            # Update sv
            # Every delay block moves one block down. The last one is deleted
            sv[hN+1:end,:] .= @view sv[1:hN*(rings-1),:]
            mul!(svu,P,prob.S.(V)) # Apply the Real Fourier Transform to V
            sv[1:hN,:].= svu # Update sv first delay ring
        end #end time loop
    end #end trajectories loop

    # Compute the mean solution across all trajectories
    @inbounds for j = 1:N
        @inbounds for i = 1:N*length(saveat)
            meanVj[i,j] = mean(@view Vj[i,j,:])
        end
    end

    # Build our solution structure output
    sol = SolveOutSto2D(Vj,meanVj,prob.Ω.x,prob.Ω.y,prob.Ω.t,saveat,N)
    
    return sol

end