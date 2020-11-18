function solveSNFE(prob,saveat=[prob.Ω.t[2], prob.Ω.t[end÷2], prob.Ω.t[end]])
    N  = prob.Ω.N
    np = prob.np
    t  = prob.Ω.t
    x  = prob.Ω.x
    y  = prob.Ω.y
    
    # Pre-allocate matrices
    Vtj = Array{Float64,3}(undef,N*length(saveat),N,np)
    meanVtj = Array{Float64,2}(undef,N*length(saveat),N)
    
    # Trajectories loop
    for p = 1:np
        j = 0 # Index for tj
        # Pre-allocate matrices
        V = Array{Float64,2}(undef,N,N)
        S = Array{Complex{Float64},2}(undef,prob.rings*N,N)
        # Pre-allocate dV..
        # Pre-allocate L..
        # Pre-allocate I (?)
        copy!(V,prob.V0) # Initial condition V(X,0) = V0
        copy!(S,prob.S)  # Initialise with the initial values of S

        # Time loop
        for i = 0:length(t)-1
            w = rand(Normal(),N,N) # Stochastic term
            t_i1=t[i+1]
            I = [prob.I(i,j,t_i1) for j in y, i in x]

            # Store solution V in t[j] instant
            if t[i+1] in saveat
                j += 1
                Vtj[(j-1)*N+1:N*j,:,p] = V
            end

            # Aqui tenho uma instabilidade de tipo. O L começa por ser complexo e depois real!!
            #usar outra matriz para armazenar a parte real de L?!
            L = prob.Krings[1:N,:].*S[1:N,:]
            for j = 2:prob.rings
                L += prob.Krings[1+N*(j-1):N*j,:].*S[1+N*(j-1):N*j,:]
            end
            # Apply the Inverse Fourier Transform
            L  = real(fftshift(ifft(ifftshift(L))))

            # Compute ....
            dV = prob.Ω.dt/prob.α*(-V+L+I) + sqrt(prob.Ω.dt)*prob.ϵ*prob.λ.*w

            # Update V
            V += dV

            # Update S
            S[N+1:end,:] = S[1:N*(prob.rings-1),:]
            S[1:N,:]     = fftshift(fft(ifftshift(prob.firingRate(V))))

        end #end time loop

    end #end trajectories loop

    # Compute the mean solution across all trajectories
    for j = 1:N
        for i = 1:N*length(saveat)
            meanVtj[i,j] = mean(Vtj[i,j,:])
        end
    end

    # Build our solution structure output
    sol = SolOutput(Vtj,meanVtj,prob.Ω.x,prob.Ω.y,t,saveat)
    
    return sol

end