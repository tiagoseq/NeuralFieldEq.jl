function probSNFE(in::Input1D)
    N  = in.N
    hN = N÷2+1 # Due the real Fourier Transform, deals with half of the dimension

    Ω = Domain{typeof(in.L),typeof(in.T),typeof(N)}(in.L,N,in.T,in.n,in.v)

    P    = plan_rfft(zeros(N),flags=FFTW.MEASURE) # Real IFFT operator for vector of dim N
    Pinv = plan_irfft(zeros(ComplexF64,hN),N,flags=FFTW.MEASURE) # Real IFFT operator for hN vectors
    #P    = plan_fft(zeros(ComplexF64,N),flags=FFTW.MEASURE) # Real IFFT operator for vector of dim N
    #Pinv = plan_ifft(zeros(ComplexF64,N),flags=FFTW.MEASURE) # Real IFFT operator for hN vectors
    # Apply rfft:  Â = P * A (==mul!(Â,P,A)). Dim of Â is N÷2+1
    # Apply irfft: A = Pinv * Â (==mul!(A,P,Â))

    @inbounds K = [in.kernel(i) for i in Ω.x] # Initialise K matrix
    mul!(K,Ω.dx,K)       # Multiply by the step size due to the discretisation of the integral operator
    Krings = peel(P,K,Ω) # The Fourier transform of K, computed on the several rings of delay

    # Initialization of past S(V) values (from t=-Tmax to t=0, where Tmax = rings*dt)
    # Each block (rings*hN) × N stores the S(V) rfft with the corresp delayed value
    sv = Vector{Complex{Float64}}(undef,Ω.rings1D*hN)
    v0 = Vector{ComplexF64}(undef,hN)
    V0 = Vector{Float64}(undef,N)
    init!(v0,sv,V0,P,in.firingRate,in.V0,Ω) # Initialise V0 in natural domain and v0 and sv in the Fourier domain

    prob = ProbOutput1D(P,Pinv,Krings,V0,v0,sv,Ω,in.α,in.extInput,in.firingRate)

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
    sv = Matrix{ComplexF64}(undef,Ω.rings2D*hN,N)
    v0 = Matrix{ComplexF64}(undef,hN,N)
    V0 = Matrix{Float64}(undef,N,N)
    init!(v0,sv,V0,P,in.firingRate,in.V0,Ω) # Initialise V0 in natural domain and v0 and sv in the Fourier domain

    prob = ProbOutput2D(P,Pinv,Krings,V0,v0,sv,Ω,in.α,in.extInput,in.firingRate)

    return prob
end