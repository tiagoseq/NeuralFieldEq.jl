function probSNFE(parameters,extIn,connectFunc,fireRate,domain)
    α, v, Vthresh, v0, ξ, ϵ, np = parameters
    L, N, tspan, n = domain 

    Ω = Domain{typeof(L),typeof(tspan)}(L,N,tspan,n)

    V0 = fill(v0,(N,N)) # initial condition V(x,y,t=0)

    # Returns our real FFT operator for NxN matrices
    P = plan_rfft(zeros(N,N))
    # Apply rfft:  Â = P * A (==mul!(Â,P,A))
    # Apply irfft: A = P \ Â (==ldiv!(A,P,Â))

    # Initialise I, K, λ matrices
    #I = [extIn(i,j) for j in Ω.y, i in Ω.x]
    @inbounds K = [connectFunc(i,j) for j in Ω.y, i in Ω.x]*Ω.dx*Ω.dy
    @inbounds λ = [exp(-((i^2+j^2)*ξ^2)/(8.0*pi)) for j in Ω.y, i in Ω.x]

    r = max(1.0,(v*Ω.dt)/Ω.dx) # compute the radius (explain what is radius)
    Krings, rings = peel(P,K,[N÷2+1,N÷2+1],r,Ω.N)

    # Initialization of past S(V) values (from t=-Tmax to t=0, where Tmax = nRings*dt)
    # Since we're working in the Fourier domain, past values are directly stored using their Fourier Transform
    fSi = Matrix{Complex{Float64}}(undef,N÷2+1,N)
    fS = repeat(mul!(fSi,P,fireRate(V0)),rings,1)

    prob = ProbOutput(P,Krings,rings,extIn,fS,fireRate,Ω,λ,α,ϵ,np,Vthresh,V0)

    # Print useful values:
    println("Spatial step: dx = dy = $(Ω.dx)")
    println("Temporal step: dt = $(Ω.dt)")
    println("Number of delay rings: $rings")

    return prob
end
