module AuxFunctions

using LinearAlgebra: ldiv!, mul!
using FFTW: fftshift
export peel, init!

function disc(centre::T,N::T,radius::P) where {T<:Signed,P<:AbstractFloat} # 1D method
    # This function constructs a disc in a vector:
    # 0 0 1 0 0 -> 0 1 0 1 0 - > 1 0 0 0 1
    D = Vector{P}(undef,N)
    @inbounds D = [sqrt((i-centre)^2)<=radius ? D[i] = 1.0 : D[i] = 0.0 for i = 1:N]
    return D
end

function disc(centre::Vector{T},N::T,radius::P) where {T<:Signed,P<:AbstractFloat} # 2D method
    # This function constructs a disc in a matrix, ex:
    # 0 0 1 0 0
    # 0 1 0 1 0
    # 0 0 1 0 0
    D = Matrix{P}(undef,N,N)
    @inbounds D = [sqrt((i-centre[1])^2+(j-centre[2])^2)<=radius ? D[i,j] = 1.0 : D[i,j] = 0.0 for i = 1:N, j = 1:N]
    return D
end


function peel(P,A::Vector{T},Ω) where {T<:AbstractFloat} # 1D method
    # This function constructs k_delay in 1D
    # Perform the Fourier Transform of the respective delay ring

    hN = Ω.N÷2+1 # Half of dim N
    B  = Vector{ComplexF64}(undef,Ω.rings1D*hN) # Each hN block corresponds to a delay ring
    Bi = Vector{ComplexF64}(undef,hN)

    # Ring loop. For each delay ring compute rfft
    @inbounds for i = 2:Ω.rings1D
        r1 = (i-1)*Ω.Δr
        r2 = i*Ω.Δr
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
    B  = Matrix{ComplexF64}(undef,hN*Ω.rings2D,Ω.N) # Each (hN X N) block corresponds to a delay ring
    Bi = Matrix{ComplexF64}(undef,hN,Ω.N)

    # Ring loop. For each delay ring compute rfft
    @inbounds for i = 2:Ω.rings2D
        r1 = (i-1)*Ω.Δr
        r2 = i*Ω.Δr
        mul!(Bi,P,fftshift((disc([hN,hN],Ω.N,r2) - disc([hN,hN],Ω.N,r1)).*A))
        B[1+hN*(i-1):hN*i,:] .= Bi
    end
    
    mul!(Bi,P,fftshift(disc([hN,hN],Ω.N,Ω.Δr).*A))
    B[1:hN,:] .= Bi
    return B
end

# Method init for 1D with V0 being a constant
function init!(v::Vector{<:Complex{<:T}},sv::Vector{<:Complex{<:T}},V0arr::Vector{T},P,S,V0::Real,Ω) where {T<:AbstractFloat}
    svu   = Vector{ComplexF64}(undef,Ω.N÷2+1)
    V0arr.= fill(V0,Ω.N) # Initial condition V(x,y,t=0)
    mul!(v,P,V0arr)      # Initial condition in the Fourier domain

    mul!(svu,P,S.(V0arr))
    sv .= repeat(svu,Ω.rings1D)
end

# Method init for 2D with V0 being a constant
function init!(v::Matrix{<:Complex{<:T}},sv::Matrix{<:Complex{<:T}},V0arr::Matrix{T},P,S,V0::Real,Ω) where {T<:AbstractFloat}
    svu   = Matrix{ComplexF64}(undef,Ω.N÷2+1,Ω.N)
    V0arr.= fill(V0,Ω.N,Ω.N) # Initial condition V(x,y,t=0)
    mul!(v,P,V0arr)          # Initial condition in the Fourier domain

    mul!(svu,P,S.(V0arr))
    sv .= repeat(svu,Ω.rings2D,1)
end

# Method init for 1D with V0 being a function
function init!(v::Vector{<:Complex{<:T}},sv::Vector{<:Complex{<:T}},V0arr::Vector{T},P,S,V0,Ω) where {T<:AbstractFloat}
    # V0 is a function
    svu   = Vector{ComplexF64}(undef,Ω.N÷2+1)
    V0arr.= [V0(i) for i in Ω.x] # Discretise V0 (Time domain)
    mul!(v,P,V0arr)              # Compute V0 in frequency domain

    mul!(svu,P,S.(V0arr))
    sv .= repeat(svu,Ω.rings1D)
end

# Method init for 2D with V0 being a function
function init!(v::Matrix{<:Complex{<:T}},sv::Matrix{<:Complex{<:T}},V0arr::Matrix{T},P,S,V0,Ω) where {T<:AbstractFloat}
    # V0 is a function
    svu   = Matrix{ComplexF64}(undef,Ω.N÷2+1,Ω.N)
    V0arr.= [V0(i,j) for j in Ω.y,i in Ω.x] # Discretise V0 (Time domain)
    mul!(v,P,V0arr)                         # Compute V0 in frequency domain

    mul!(svu,P,S.(V0arr))
    sv .= repeat(svu,Ω.rings2D,1)
end

end #end module
