@doc raw"""
    probSNFE(input)

Set the NFE problem by discretising the domain, initialising `sv` and `k_delay`
and building the delay rings structure.

Receive as input the parameters and functions wrapped in structures `Input1D` and `Input2D`.
According to the space dimension the inputs are wrapped in the corresponding structure.
The arguments of `Input1D` and `Input2D` are:
# Arguments
- `α::AbstractFloat`
- `v::AbstractFloat`
- `V0::fV0`: Can be a number or a function
- `L::Number`
- `N::Integer`
- `T::AbstractFloat`
- `n::Integer`
- `I::fI`: External input function
- `K::fK`: Connectivity function
- `S::fS`: Firing rate function

Return a structure to be used in solveSNFE and prints useful information about the problem.

# Examples
```julia-repl
julia> # 2D example
julia> I(x,y,t) = 0                      # External inut
julia> # In 2D even if I doesn't depend on x,y and t, these three arguments are mandatory
julia> K(x,y)   = exp(-x^2+y^2)          # Connectivity function
julia> S(V)     = convert(Float64,V>0.0) # Firing rate. Heavyside function H(V)

julia> α  = 1.0   # Constant decay      (must be float)
julia> v  = 20.0  # Finite axonal speed (must be float)
julia> V0 = 0.0   # Initial condition (can be a constant or a function)
julia> L  = 100   # Domain length     (can be a integer or float)
julia> N  = 512   # Number of nodes to discretise space (must be integer)
julia> T  = 20.0  # Time span (must be float)
julia> n  = 200   # Number of nodes to discretise time  (must be integer)

julia> input = Input2D(α,v,V0,L,N,T,n,I,K,S);
julia> prob  = probSNFE(input)
├─ Domain:       Ω × [0,T] = [-50.0,49.8046875]^2 × [0,20.0]
├─ Spatial step: dx   = 0.1953125
├─ Time step:    dt   = 0.1
├─ Velocity:     v    = 20.0
├─ Delay rings:  umax = 36
```

For the 1D scenario the procedure is the same, with exception that the functions
`I` and `K` must be declared as `I(x,t)` and `K(x)` using `Input1D`.
"""
function probSNFE(in::Input1D)
    N  = in.N
    hN = N÷2+1 # Due the real Fourier Transform, deals with half of the dimension

    # Contruct our discretised domain
    Ω = Domain{typeof(in.L),typeof(in.T),typeof(N)}(in.L,N,in.T,in.n,in.v)

    P    = plan_rfft(zeros(N),flags=FFTW.MEASURE) # Real IFFT operator for vector of dim N
    Pinv = plan_irfft(zeros(ComplexF64,hN),N,flags=FFTW.MEASURE) # Real IFFT operator for hN vectors
    # Apply rfft:  Â = P * A (==mul!(Â,P,A)). Dim of Â is N÷2+1
    # Apply irfft: A = Pinv * Â. Do not use mul! in the irfft case

    @inbounds K = [in.kernel(i) for i in Ω.x] # Initialise K matrix
    mul!(K,Ω.dx,K)       # Multiply K by the step size due to the discretisation of the integral operator
    krings = peel(P,K,Ω) # The rfft of K computed on the several rings of delay (k_delay)

    # Initialization of past S(V) values (from t=-Tmax to t=0, where Tmax = rings*dt)
    # Each delay ring block (hN × N) stores the S(V) rfft with the corresp delayed value
    sv = Vector{ComplexF64}(undef,Ω.rings1D*hN)
    v0 = Vector{ComplexF64}(undef,hN)
    V0 = Vector{Float64}(undef,N)
    # Initialise V0 in time domain and v0 and sv in the Fourier domain
    init!(v0,sv,V0,P,in.firingRate,in.V0,Ω)

    prob = ProbOutput1D(P,Pinv,krings,V0,v0,sv,Ω,in.α,in.extInput,in.firingRate)

return prob
end

function probSNFE(in::Input2D)
    N  = in.N
    hN = N÷2+1 # Due the real Fourier Transform, deals with half of the dimension

    # Contruct our discretised domain
    Ω = Domain{typeof(in.L),typeof(in.T),typeof(N)}(in.L,N,in.T,in.n,in.v)

    P    = plan_rfft(zeros(N,N),flags=FFTW.MEASURE) # Real FFT operator for NxN matrices
    Pinv = plan_irfft(zeros(ComplexF64,hN,N),N,flags=FFTW.MEASURE) # Real IFFT operator for hNxN matrices
    # Apply rfft:  Â = P * A (==mul!(Â,P,A)). Dim of Â is ((N÷2+1) × N)
    # Apply irfft: A = Pinv * Â. Do not use mul! in the irfft case

    @inbounds K = [in.kernel(i,j) for j in Ω.y, i in Ω.x] # Initialise K matrix
    mul!(K,Ω.dx*Ω.dy,K)  # Multiply K by the step size due to the discretisation of the integral operator
    krings = peel(P,K,Ω) # The rfft of K computed on the several rings of delay (k_delay)

    # Initialization of past S(V) values (from t=-Tmax to t=0, where Tmax = rings*dt)
    # Each delay ring block (rings*hN) × N stores the S(V) rfft with the corresp delayed value
    sv = Matrix{ComplexF64}(undef,Ω.rings2D*hN,N)
    v0 = Matrix{ComplexF64}(undef,hN,N)
    V0 = Matrix{Float64}(undef,N,N)
    # Initialise V0 in time domain and v0 and sv in the Fourier domain
    init!(v0,sv,V0,P,in.firingRate,in.V0,Ω)

    prob = ProbOutput2D(P,Pinv,krings,V0,v0,sv,Ω,in.α,in.extInput,in.firingRate)

    return prob
end