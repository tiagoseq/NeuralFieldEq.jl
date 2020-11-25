function solveSNFE(prob,saveat=[prob.Ω.t[2], prob.Ω.t[end÷2], prob.Ω.t[end]])
    P  = prob.P
    N  = prob.Ω.N
    np = prob.np
    t  = prob.Ω.t
    x  = prob.Ω.x
    y  = prob.Ω.y

    hN = (N÷2)+1 # Half of dim N (due to the output of rfft)
    
    # Pre-allocate the solution matrices
    Vtj = Array{Float64,3}(undef,N*length(saveat),N,np)
    meanVtj = Array{Float64,2}(undef,N*length(saveat),N)
    
    # Trajectories loop
    @inbounds for p = 1:np
        j = 0 # Index for tj

        # Pre-allocate matrices
        V   = Array{Float64,2}(undef,N,N)
        dv  = similar(V)
        L   = similar(V)
        fS  = Array{Complex{Float64},2}(undef,prob.rings*hN,N) # Firing Rate in Fourier space
        fSi = Array{Complex{Float64},2}(undef,hN,N) # Firing Rate in Fourier space
        fL  = similar(fSi) # kernel by firing rate product in Fourier space
        # Pre-allocate I (?)

        copy!(V,prob.V0) # Initial condition V(X,0) = V0
        copy!(fS,prob.S)  # Initialise with the initial values of S

        # Time loop
        @inbounds for i = 0:length(t)-1
            w = rand(Normal(),N,N) # Stochastic term
            t_i1=t[i+1]
            I = [prob.I(i,j,t_i1) for j in y, i in x]

            # Store solution V in t[j] instant
            if t[i+1] in saveat
                j += 1
                Vtj[(j-1)*N+1:N*j,:,p] = V
            end

            fL = prob.Krings[1:hN,:].*fS[1:hN,:]
            @inbounds for j = 2:prob.rings
                fL += prob.Krings[1+hN*(j-1):hN*j,:].*fS[1+hN*(j-1):hN*j,:]
            end

            # Apply the Inverse Real Fourier Transform
            ldiv!(L,P,fL)

            # Compute dV - explain what is dV
            dV = (prob.Ω.dt/prob.α)*(-V+L+I) + sqrt(prob.Ω.dt)*prob.ϵ*prob.λ.*w

            V += dV # Update V

            # Update S
            fS[hN+1:end,:] = fS[1:hN*(prob.rings-1),:]
            
            mul!(fSi,P,prob.firingRate(V))
            fS[1:hN,:] .= fSi

        end #end time loop

    end #end trajectories loop

    # Compute the mean solution across all trajectories
    @inbounds for j = 1:N
        @inbounds for i = 1:N*length(saveat)
            meanVtj[i,j] = mean(Vtj[i,j,:])
        end
    end

    # Build our solution structure output
    sol = SolOutput(Vtj,meanVtj,prob.Ω.x,prob.Ω.y,t,saveat)

    println("Solution saved at: $(saveat) seconds")
    
    return sol

end