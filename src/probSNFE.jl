function probSNFE(in::Input1D)
    N  = in.N
    hN = N÷2+1 # Due the real Fourier Transform, deals with half of the dimension

    Ω = Domain{typeof(in.L),typeof(in.T),typeof(N)}(in.L,N,in.T,in.n,in.v)

    P    = plan_rfft(zeros(N)) # Real IFFT operator for vector of dim N
    Pinv = plan_irfft(zeros(ComplexF64,hN),N,flags=FFTW.MEASURE) # Real IFFT operator for hN vectors
    # Apply rfft:  Â = P * A (==mul!(Â,P,A)). Dim of Â is N÷2+1
    # Apply irfft: A = Pinv * Â (==mul!(A,P,Â))

    @inbounds K = [in.kernel(i) for i in Ω.x] # Initialise K matrix
    mul!(K,Ω.dx,K)       # Multiply by the step size due to the discretisation of the integral operator
    Krings = peel(P,K,Ω) # The Fourier transform of K, computed on the several rings of delay

    # Initialization of past S(V) values (from t=-Tmax to t=0, where Tmax = rings*dt)
    # Each block (rings*hN) × N stores the S(V) rfft with the corresp delayed value
    s  = Vector{Complex{Float64}}(undef,Ω.rings*hN)
    v0 = Vector{ComplexF64}(undef,hN)
    init!(v0,s,P,in.firingRate,in.V0,Ω) # Initialise v0 and s in the Fourier domain

    prob = ProbOutput1D(P,Krings,v0,s,Ω,in.α,in.extInput,in.firingRate)

    # Print useful values:
    println("Spatial step: dx = $(Ω.dx)")
    println("Temporal step: dt = $(Ω.dt)")
    println("Ω × [0;T] = [$(Ω.x[1]),$(Ω.x[end])] × [0, $(Ω.t[end])]")
    println("Velocity: $(in.v)")
    println("Number of delay rings: $(Ω.rings)")

    return prob
end

function probSNFE(in::Input2D)
    N  = in.N
    hN = N÷2+1 # Due the real Fourier Transform, deals with half of the dimension

    Ω = Domain{typeof(in.L),typeof(in.T),typeof(N)}(in.L,N,in.T,in.n,in.v)

    P    = plan_rfft(zeros(N,N),flags=FFTW.MEASURE) # Real FFT operator for NxN matrices
    Pinv = plan_irfft(zeros(ComplexF64,hN,N),N,flags=FFTW.MEASURE) # Real IFFT operator for hNxN matrices
    # Apply rfft:  Â = P * A (==mul!(Â,P,A)). Dim of Â is ((N÷2+1) × N)
    # Apply irfft: A = Pinv * Â (==mul!(A,P,Â))

    @inbounds K = [in.kernel(i,j) for j in Ω.y, i in Ω.x] # Initialise K matrix
    mul!(K,Ω.dx*Ω.dy,K)  # Multiply by the step size due to the discretisation of the integral operator
    Krings = peel(P,K,Ω) # The Fourier transform of K computed on the several rings of delay

    # Initialization of past S(V) values (from t=-Tmax to t=0, where Tmax = rings*dt)
    # Each block (rings*hN) × N stores the S(V) rfft with the correspondent delayed value
    s  = Matrix{ComplexF64}(undef,Ω.rings*hN,N)
    v0 = Matrix{ComplexF64}(undef,hN,N)
    init!(v0,s,P,in.firingRate,in.V0,Ω) # Initialise v0 and s in the Fourier domain

    prob = ProbOutput2D(P,Pinv,Krings,v0,s,Ω,in.α,in.extInput,in.firingRate)

    # Print useful values:
    println("Spatial step: dx = dy = $(Ω.dx)")
    println("Temporal step: dt = $(Ω.dt)")
    println("Ω × [0;T] = [$(Ω.x[1]),$(Ω.x[end])]^2 × [0, $(Ω.t[end])]")
    println("Velocity: $(in.v)")
    println("Number of delay rings: $(Ω.rings)")

    return prob
end