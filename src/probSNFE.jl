function probSNFE(in::Input1D)
    v  = in.v
    v0 = in.v0
    L  = in.L
    N  = in.N
    T  = in.T
    n  = in.n
    hN = N÷2+1 # Due the real Fourier Transform, deals with half of the dimension

    Ω = Domain{typeof(L),typeof(T)}(L,N,T,n)

    V0 = fill(v0,N) # initial condition V(x,t=0)

    # Returns our real FFT operator for vector of dim N
    P = plan_rfft(zeros(N))
    # Apply rfft:  Â = P * A (==mul!(Â,P,A)). Dim of Â is N÷2+1
    # Apply irfft: A = P \ Â (==ldiv!(A,P,Â))

    @inbounds K = [in.kernel(i) for i in Ω.x]*Ω.dx # Initialise K matrix

    r = max(1.0,(v*Ω.dt)/Ω.dx) # compute the radius (explain what is radius)
    Krings, rings = peel(P,K,hN,r,Ω.N)

    # Initialization of past S(V) values (from t=-Tmax to t=0, where Tmax = nRings*dt)
    # Since we're working in the Fourier domain, past values are directly stored using their Fourier Transform
    si = Vector{Complex{Float64}}(undef,hN)
    s  = Vector{Complex{Float64}}(undef,rings*hN) # Firing Rate in Fourier space
    s  = repeat(mul!(si,P,in.firingRate(V0)),rings,1)

    prob = ProbOutput1D(P,Krings,rings,s,V0,Ω,in)

    # Print useful values:
    println("Spatial step: dx = $(Ω.dx)")
    println("Temporal step: dt = $(Ω.dt)")
    println("Ω × [0;T] = [$(Ω.x[1]),$(Ω.x[end])] × [0, $(Ω.t[end])]")
    println("Number of delay rings: $rings")

    return prob
end

function probSNFE(in::Input2D)
    v  = in.v
    v0 = in.V0
    L  = in.L
    N  = in.N
    T  = in.T
    n  = in.n
    hN = N÷2+1 # Due the real Fourier Transform, deals with half of the dimension

    Ω = Domain{typeof(L),typeof(T)}(L,N,T,n)

    V0 = fill(v0,(N,N)) # initial condition V(x,y,t=0)

    # Returns our real FFT operator for NxN matrices
    P = plan_rfft(zeros(N,N))
    #P = plan_fft(zeros(N,N))
    # Apply rfft:  Â = P * A (==mul!(Â,P,A)). Dim of Â is ((N÷2+1) × N)
    # Apply irfft: A = P \ Â (==ldiv!(A,P,Â))

    @inbounds K = [in.kernel(i,j) for j in Ω.y, i in Ω.x]*Ω.dx*Ω.dy # Initialise K matrix

    r = max(1.0,(v*Ω.dt)/Ω.dx) # compute the radius (explain what is radius)
    Krings, rings = peel(P,K,[N÷4+1,hN],r,Ω.N)

    # Initialization of past S(V) values (from t=-Tmax to t=0, where Tmax = nRings*dt)
    # Since we're working in the Fourier domain, past values are directly stored using their Fourier Transform
    si = Matrix{Complex{Float64}}(undef,hN,N)
    s  = Matrix{Complex{Float64}}(undef,rings*hN,N) # Firing Rate in Fourier space
    s  = repeat(mul!(si,P,in.firingRate(V0)),rings,1)

    prob = ProbOutput2D(P,Krings,rings,s,V0,Ω,in)

    # Print useful values:
    println("Spatial step: dx = dy = $(Ω.dx)")
    println("Temporal step: dt = $(Ω.dt)")
    println("Ω × [0;T] = [$(Ω.x[1]),$(Ω.x[end])]^2 × [0, $(Ω.t[end])]")
    println("Number of delay rings: $rings")

    return prob
end