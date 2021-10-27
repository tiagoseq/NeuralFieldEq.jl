@doc raw"""
    solveNFE(problem,saveat)

# Arguments:
- `problem :: ProbOutput1D`: Output of probNFE using Input1D structure
- `problem :: ProbOutput2D`: Output of probNFE using Input2D structure
- `saveat  :: AbstractVector`: Vector containing the instants where the solution is saved at

Returns a structure containing the solution to the NFE saved at saveat instants.
The structure has `x`, `y` and `saveat` as fields to help plotting the solution.

# Examples
```julia-repl
julia> # 2D example
julia> saveat = [5.0,20.0]
julia> sol    = solveNFE(prob,saveat) # Deterministic solution saved at t=5 and t=20
julia> sol(20.0) # Solution at t=20
julia> x = sol.x # Spatial vector in x direction 
julia> y = sol.y # Spatial vector in y direction
julia> using Plots
julia> plot(x,y,sol(20.0),st=:surface,title="Solution at t=20") # Surface plot of sol(20.0)
```
If the neural field is in 1D/2D, the solution is a vector/matrix.

-------------------------------------------------------------------------------------------
    solveNFE(problem,saveat,ϵ,np,ξ=0.1)

Solve stochastic version of an NFE for np trajectories, noise level ϵ and correlation ξ.

# Arguments:
- `problem :: ProbOutput1D`: Output of probNFE with Input1D 
- `problem :: ProbOutput2D`: Output of probNFE with Input2D
- `saveat  :: AbstractVector`: Vector containing the instants where the solution is saved at
- `ϵ       :: Number`: Level of additive noise
- `np      :: Integer`: Number of simulations
- `ξ       :: AbstractFloat=0.1`: Spatial correlation parameter

Returns a structure containing the mean solution and the trajectories of the
NFE saved at saveat instants.

# Examples
```julia-repl
julia> # Stochastic solution saved at t=5,20, with noise level 0.01 simulated 100 times.
julia> sol_sto = solveNFE(prob,saveat,0.01,100)
julia> sol_sto(20.0,4) # 4th path at t=20
julia> sol_sto(5.0)    # Mean solution at t=20
```
"""
function solveNFE(prob::ProbOutput1D,saveat) # Deterministic method for 1D
    P      = prob.Plan
    Pinv   = prob.PlanInv
    krings = prob.krings
    α      = prob.α
    rings  = prob.Ω.rings1D
    N      = prob.Ω.N
    t      = prob.Ω.t
    x      = prob.Ω.x
    dt     = prob.Ω.dt
    V      = prob.V0 # Initial condition in natural space (needed if saveat contains t[1])
    v      = prob.v0 # Initial condition in Fourier space
    sv     = prob.sv # Firing Rate in Fourier space (all delays)
    hN     = N÷2+1   # Half of dim N (due to the output of rfft)

    # Pre-allocate vectors and matrices
    Vj  = Vector{Float64}(undef,N*length(saveat)) # Store the solution at saveat
    svu = Vector{ComplexF64}(undef,hN) # Firing Rate in Fourier space (delay u)
    I   = similar(V)
    î   = similar(svu) # rfft of I
    a   = similar(svu) # rfft of the integral operator A(X,t)
    j   = 1 # Index for saveat

    @inbounds for i = 1:length(t) # Time loop
        ti = t[i]
        # Store solution V at saveat instants
        if ti in saveat
            Vj[(j-1)*N+1:N*j] .= V
            j += 1
        end

        I.= [prob.I(i,ti) for i in x] # Discretise I
        mul!(î,P,I)                   # Compute rfft of I

        # Compute the integral operator's at frequency domain at t_i
        @. a = @views krings[1:hN]*sv[1:hN] # Init a
        @inbounds for u = 2:rings # Rings loop
            @. a += @views krings[1+hN*(u-1):hN*u]*sv[1+hN*(u-1):hN*u]
        end

        # Euler explicit scheme: v_i+1 = v_i + (Δt/α)*(î_i - v_i + a_i)
        @. v = v + (dt/α)*(î - v + a)

        V .= Pinv * v # Perform inverse Fourier transform to compute V in natural domain

        # Update sv
        # Every delay block moves one block down. The last one is deleted
        sv[hN+1:end] .= @view sv[1:hN*(rings-1)] 
        mul!(svu,P,prob.S.(V)) # Apply the Real Fourier Transform to V
        sv[1:hN] .= svu        # Update sv first delay ring
    end #end time loop

    # Build our solution structure output
    sol = SolveOutDet1D(Vj,x,t,saveat,N)
    
    return sol

end

# Stochastic method for 1D
function solveNFE(prob::ProbOutput1D,saveat,ϵ,np,ξ=0.1)
    P      = prob.Plan
    Pinv   = prob.PlanInv
    krings = prob.krings
    α      = prob.α
    rings  = prob.Ω.rings1D
    N      = prob.Ω.N
    t      = prob.Ω.t
    x      = prob.Ω.x
    dt     = prob.Ω.dt
    hN     = N÷2+1   # Half of dim N (due to the output of rfft)
    
    # Pre-allocate matrices
    meanVj = Vector{Float64}(undef,N*length(saveat)) # Store the mean solution at saveat
    Vj  = Matrix{Float64}(undef,N*length(saveat),np) # Store the solution at saveat at path p
    V   = Vector{Float64}(undef,N)
    sv  = Vector{ComplexF64}(undef,hN*rings) # Firing Rate in Fourier space (all delays)
    svu = Vector{ComplexF64}(undef,hN) # Firing Rate in Fourier space (delay u)
    I   = similar(V)
    v   = similar(svu) 
    î   = similar(svu) # rfft of I
    a   = similar(svu) # Fourier coefficients of the integral operator A(X,t)
    λ   = [exp(-(i^2*ξ^2)/(8*pi)) for i = 1:hN] # Correlation matrix
    
    # Pre-compute constant term of the stochastic part
    # fft gives the fourier coeff mult by N. Need to scale up noise, to have same magnitude
    # √dt due to the explicit Euler-Maruyama
    # Shift λ due to the also shifted output of fft
    c = N*sqrt(dt)*ϵ*fftshift(λ) # ϵ is the level of additive noise
    
    # Trajectories loop
    @inbounds for p = 1:np
        copy!(V,prob.V0)  # Initial condition V(X,0) = V0 (Natural domain)
        copy!(v,prob.v0)  # Initial condition v(x,0) = v0 (Fourier domain)
        copy!(sv,prob.sv) # Initialise with the initial values of sv
        
        j = 1 # Index for tj
        @inbounds for i = 1:length(t) # Time loop
            # w ~ N(0,1) + i*N(0,1)  stochastic term
            w = randn(Float64,hN) + im*randn(Float64,hN)
            ti = t[i]
            # Store solution V in saveat instants
            if ti in saveat
                Vj[(j-1)*N+1:N*j,p] .= V
                j += 1
            end
    
            I.= [prob.I(i,ti) for i in x] # Discretise I
            mul!(î,P,I)                   # rfft I
    
            # Compute the integral operator's at frequency domain at t_i
            @. a = @views krings[1:hN]*sv[1:hN] # Init a
            @inbounds for u = 2:rings # Ring loop
                @. a += @views krings[1+hN*(u-1):hN*u]*sv[1+hN*(u-1):hN*u]
            end

            # Euler-Maruyama explicit scheme: v_i+1 = v_i + (Δt/α)*(î_i - v_i + a_i) + √(Δt)*ϵ*λ*w
            @. v = v + (dt/α)*(î - v + a) + c*w

            V .= Pinv * v # Perform inverse Fourier transform to compute V in natural domain

            # Update sv
            # Every delay block moves one block down. The last one is deleted
            sv[hN+1:end] .= @view sv[1:hN*(rings-1)]
            mul!(svu,P,prob.S.(V)) # Apply the Real Fourier Transform to V
            sv[1:hN,:].= svu       # Update sv first delay ring
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