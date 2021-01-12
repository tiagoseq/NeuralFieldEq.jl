# Deterministic method for 1D
function solveSNFE(prob::ProbOutput1D,saveat)
    P  = prob.P
    α  = prob.in.α
    N  = prob.Ω.N
    t  = prob.Ω.t
    x  = prob.Ω.x
    dt = prob.Ω.dt
    s  = prob.s # Firing rate history in Fourier Space (vector)
    hN = (N÷2)+1 # Half of dim N (due to the output of rfft)
    
    # Pre-allocate vectors and matrices
    Vtj = Vector{Float64}(undef,N*length(saveat))
    V   = Vector{Float64}(undef,N)
    dV  = similar(V)
    L   = similar(V)
    I   = similar(V)
    si  = Vector{Complex{Float64}}(undef,hN) # Firing Rate in Fourier space (delay i)
    l   = similar(si) # kernel by firing rate product in Fourier space
    copy!(V,prob.V0) # Initial condition V(X,0) = V0
    copy!(s,prob.s)  # Initialise with the initial values of S
    j = 0 # Index for tj

    # Time loop
    @inbounds for i = 0:length(t)-1
        t_i1=t[i+1]
        I = [prob.in.extInput(i,t_i1) for i in x]

        # Store solution V in t[j] instant
        if t[i+1] in saveat
            j += 1
            Vtj[(j-1)*N+1:N*j] .= V
        end

        @. l = prob.Krings[1:hN,:]*s[1:hN,:]
        @inbounds for j = 2:prob.rings
            @. l += prob.Krings[1+hN*(j-1):hN*j,:]*s[1+hN*(j-1):hN*j,:]
        end

        # Apply the Real Inverse Fourier Transform
        ldiv!(L,P,l)

        @. dV = (dt/α)*(-V+L+I) # Compute dV -> explain what is dV!!
        @. V += dV # Update V in time

        # Update s
        s[hN+1:end,:] .= s[1:hN*(prob.rings-1),:]
        mul!(si,P,prob.in.firingRate(V)) # Apply Real Fourier Transform
        s[1:hN,:] .= si

    end #end time loop

    # Build our solution structure output
    sol = SolveOutDet1D(Vtj,x,y,t,saveat)

    println("Solution saved at time instants: $(saveat)")
    
    return sol

end

# Stochastic method for 1D
function solveSNFE(prob::ProbOutput1D,saveat,(ϵ,ξ,np))
    P  = prob.P
    α  = prob.in.α
    N  = prob.Ω.N
    t  = prob.Ω.t
    x  = prob.Ω.x
    y  = prob.Ω.y
    dt = prob.Ω.dt
    hN = (N÷2)+1 # Half of dim N (due to the output of rfft)
    
    # Pre-allocate matrices
    Vtj = Matrix{Float64}(undef,N*length(saveat),np)
    meanVtj = Matrix{Float64}(undef,N*length(saveat),N)
    V  = Matrix{Float64}(undef,N,N)
    dV = similar(V)
    L  = similar(V)
    I  = similar(V)
    si = Matrix{Complex{Float64}}(undef,hN,N) # Firing Rate in Fourier space (delay i)
    l  = similar(si) # kernel ring by firing rate product in Fourier space
    λ  = [exp(-((i^2+j^2)*ξ^2)/(8*pi)) for j in y, i in x] # Correlation matrix
    
    # Trajectories loop
    @inbounds for p = 1:np
        copy!(V,prob.V0) # Initial condition V(X,0) = V0
        copy!(s,prob.s)  # Initialise with the initial values of S
        j = 0 # Index for tj

        # Time loop
        @inbounds for i = 0:length(t)-1
            w = rand(Normal(),N,N) # Stochastic term
            t_i1=t[i+1]
            I = [prob.in.extInput(i,j,t_i1) for j in y, i in x]

            # Store solution V in t[j] instant
            if t[i+1] in saveat
                j += 1
                Vtj[(j-1)*N+1:N*j,p] .= V
            end

            @. l = prob.Krings[1:hN,:]*s[1:hN,:]
            @inbounds for j = 2:prob.rings
                @. l += prob.Krings[1+hN*(j-1):hN*j,:]*s[1+hN*(j-1):hN*j,:]
            end

            # Apply the Real Inverse Fourier Transform
            ldiv!(L,P,l)

            @. dV = (dt/α)*(-V+L+I) + sqrt(dt)*ϵ*λ*w # Compute dV -> explain what is dV!!
            @. V += dV # Update V in time

            # Update s
            s[hN+1:end,:] .= s[1:hN*(prob.rings-1),:]
            mul!(si,P,prob.in.firingRate(V))
            s[1:hN,:] .= si

        end #end time loop

    end #end trajectories loop

    # Compute the mean solution across all trajectories
    @inbounds for i = 1:N*length(saveat)
        meanVtj[i] = mean(Vtj[i,:])
    end

    # Build our solution structure output
    sol = SolveOutSto1D(Vtj,meanVtj,x,y,t,saveat)

    println("Solution saved at time instants: $(saveat)")
    
    
    return sol

end