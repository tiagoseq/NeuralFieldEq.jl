@doc raw"""
This module provides auxiliary functions used internally in our method.
Is divided into two sections. Regarding to the computation of k_delay,
functions: `disc` and `peel`. Regarding to the initialisation of V0 and S(V0),
function: `init` (with several methods).

### `disc(centre,N,radius)`:
Constructs a compact disc in vectors and matrices of size `N` and `N x (N÷2+1)`, respectively.
Examples:
- Vector with `centre=3`, size: `N=5`, `radius=0` -> 0 0 1 0 0
- Vector with `centre=3`, size: `N=5`, `radius=1` -> 0 1 1 1 0
- Vector with `centre=3`, size: `N=5`, `radius=2` -> 1 1 1 1 1
- Matrix with `centre=([2,2])`, size: `N=3`, `radius=1` -> [0 1 0; 1 1 1; 0 1 0]
To call this function: `NeuralFieldEq.AuxFunctions.disc`.

### `peel(P,A,Ω)`:
Constructs the matrix k_delay (in Fourier domain).
Based on function disc computes the fft for every concentric ring present in Ω.
The fft of the delay ring `j` is stored at the jth row (in 1D) and jth+N rows (in 2D).
To call this function: `NeuralFieldEq.AuxFunctions.peel`.

### `init!(v,sv,V0arr,P,S,V0,Ω)`:
Initialises in Fourier domain the initial condition `V0` and the image of `S(V0)`.
Returns a vector (size N) According to the dimensionality of the problem.
To call this function: `NeuralFieldEq.AuxFunctions.init!`.
"""
module AuxFunctions

using LinearAlgebra: mul!
using FFTW: fftshift
export peel, init!

function disc(centre::T,N::T,radius::P) where {T<:Signed,P<:AbstractFloat} # 1D method
    D = Vector{P}(undef,N)
    @inbounds D = [sqrt((i-centre)^2)<=radius ? D[i] = 1.0 : D[i] = 0.0 for i = 1:N]
    return D
end

function disc(centre::Vector{T},N::T,radius::P) where {T<:Signed,P<:AbstractFloat} # 2D method
    D = Matrix{P}(undef,N,N)
    @inbounds D = [sqrt((i-centre[1])^2+(j-centre[2])^2)<=radius ? D[i,j] = 1.0 : D[i,j] = 0.0 for i = 1:N, j = 1:N]
    return D
end


function peel(P,A::Vector{T},Ω) where {T<:AbstractFloat} # 1D method
    # This function constructs k_delay in 1D
    # Perform the Fourier Transform of the respective delay ring

    hN = Ω.N÷2+1 # Half of dim N
    B  = Vector{Complex{T}}(undef,Ω.rings1D*hN) # Each hN block corresponds to a delay ring
    Bi = Vector{Complex{T}}(undef,hN)

    # Ring loop. For each delay ring compute rfft
    @inbounds for i = 2:Ω.rings1D
        r1 = (i-1)*Ω.Δr
        r2 = i*Ω.Δr
        # fftshift is needed because the integral op A is a circulary convolution
        mul!(Bi,P,fftshift((disc(hN,Ω.N,r2) - disc(hN,Ω.N,r1)).*A))
        B[1+hN*(i-1):hN*i] .= Bi
    end

    mul!(Bi,P,fftshift(disc(hN,Ω.N,Ω.Δr).*A))
    B[1:hN] .= Bi
    return B
end

function peel(P,A::Matrix{T},Ω) where {T<:AbstractFloat} # 2D method
    # This function constructs k_delay in 2D
    # Perform the Fourier Transform of the respective delay ring

    hN = Ω.N÷2+1 # Half of dim N
    B  = Matrix{Complex{T}}(undef,hN*Ω.rings2D,Ω.N) # Each (hN X N) block corresponds to a delay ring
    Bi = Matrix{Complex{T}}(undef,hN,Ω.N)

    # Ring loop. For each delay ring compute rfft
    @inbounds for i = 2:Ω.rings2D
        r1 = (i-1)*Ω.Δr
        r2 = i*Ω.Δr
        # fftshift is needed because the integral op A is a circulary convolution
        mul!(Bi,P,fftshift((disc([hN,hN],Ω.N,r2) - disc([hN,hN],Ω.N,r1)).*A))
        B[1+hN*(i-1):hN*i,:] .= Bi
    end
    
    mul!(Bi,P,fftshift(disc([hN,hN],Ω.N,Ω.Δr).*A))
    B[1:hN,:] .= Bi
    return B
end

# Method init for 1D with V0 being a constant
function init!(v::Vector{<:Complex{<:T}},sv::Vector{<:Complex{<:T}},V0arr::Vector{T},P,S,V0::Real,Ω) where {T<:AbstractFloat}
    svu   = Vector{Complex{T}}(undef,Ω.N÷2+1)
    V0arr.= fill(V0,Ω.N) # Initial condition V(x,y,t=0)
    mul!(v,P,V0arr)      # Initial condition in the Fourier domain

    mul!(svu,P,S.(V0arr))
    sv .= repeat(svu,Ω.rings1D)
end

# Method init for 2D with V0 being a constant
function init!(v::Matrix{<:Complex{<:T}},sv::Matrix{<:Complex{<:T}},V0arr::Matrix{T},P,S,V0::Real,Ω) where {T<:AbstractFloat}
    svu   = Matrix{Complex{T}}(undef,Ω.N÷2+1,Ω.N)
    V0arr.= fill(V0,Ω.N,Ω.N) # Initial condition V(x,y,t=0)
    mul!(v,P,V0arr)          # Initial condition in the Fourier domain

    mul!(svu,P,S.(V0arr))
    sv .= repeat(svu,Ω.rings2D,1)
end

# Method init for 1D with V0 being a function
function init!(v::Vector{<:Complex{<:T}},sv::Vector{<:Complex{<:T}},V0arr::Vector{T},P,S,V0,Ω) where {T<:AbstractFloat}
    # V0 is a function
    svu   = Vector{Complex{T}}(undef,Ω.N÷2+1)
    V0arr.= [V0(i) for i in Ω.x] # Discretise V0 (Time domain)
    mul!(v,P,V0arr)              # Compute V0 in frequency domain

    mul!(svu,P,S.(V0arr))
    sv .= repeat(svu,Ω.rings1D)
end

# Method init for 2D with V0 being a function
function init!(v::Matrix{<:Complex{<:T}},sv::Matrix{<:Complex{<:T}},V0arr::Matrix{T},P,S,V0,Ω) where {T<:AbstractFloat}
    # V0 is a function
    svu   = Matrix{Complex{T}}(undef,Ω.N÷2+1,Ω.N)
    V0arr.= [V0(i,j) for j in Ω.y,i in Ω.x] # Discretise V0 (Time domain)
    mul!(v,P,V0arr)                         # Compute V0 in frequency domain

    mul!(svu,P,S.(V0arr))
    sv .= repeat(svu,Ω.rings2D,1)
end

end #end module
