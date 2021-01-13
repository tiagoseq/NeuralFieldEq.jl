module AuxFunctions

using LinearAlgebra: ldiv!, mul!
export peel

function disc(N,centre::Int64,radius) # disc method for 1D
    # This function constructs a disc in a vector:
    # 0 0 1 0 0 -> 0 1 0 1 0 - > 1 0 0 0 1
    D = Vector{Float64}(undef,N)
    @inbounds D = [sqrt((i-centre)^2)<=radius ? D[i] = 1.0 : D[i] = 0.0 for i = 1:N]

    return D
end

function disc(N,centre::Vector{Int64},radius) # disc method for 2D
    # This function constructs a disc in a matrix:
    # 0 0 1 0 0
    # 0 1 0 1 0
    # 0 0 1 0 0
    D = Matrix{Float64}(undef,N,N)
    @inbounds D = [sqrt((i-centre[1])^2+(j-centre[2])^2)<=radius ? D[i,j] = 1.0 : D[i,j] = 0.0 for i = 1:N, j = 1:N]

    return D
end


function peel(P,Z::Vector{T},centre::M,r,N) where {T<:AbstractFloat,M<:Signed} # peel method for 1D
    # Compute the maximum diameter to get number of rings
    rx    = max(N-centre-1,centre-1)
    R     = abs(rx)
    rings = 1+floor(Int64, R/r) # Number of rings
    hN    = N÷2+1 # Half of dim N (due to the output of rfft)

    L  = Vector{Complex{Float64}}(undef,rings*hN)
    Li = Vector{Complex{Float64}}(undef,hN)

    @inbounds for i = 1:rings
        r1 = (i-1)*r
        r2 = i*r
        if i > 1
            mul!(Li,P,(disc(N,centre,r2) - disc(N,centre,r1)).*Z)
            L[1+hN*(i-1):hN*i] .= Li            
        else
            mul!(Li,P,Z)
            L[1:hN] .= Li
        end
    end

    return L, rings
end

function peel(P,Z::Matrix{T},centre::Vector{M},r,N) where {T<:AbstractFloat,M<:Signed} # peel method for 2D
    # Compute the maximum diameter to get number of rings
    rx    = max(N-centre[1],centre[1])
    ry    = max(N-centre[2],centre[2])
    R     = sqrt(rx^2+ry^2)
    rings = 1+floor(Int64, R/r) # Number of delay rings
    hN    = N÷2+1 # Half of dim N (due to the output of rfft)

    L  = Matrix{Complex{Float64}}(undef,hN*rings,N)
    Li = Matrix{Complex{Float64}}(undef,hN,N)
    @inbounds for i = 1:rings
        r1 = (i-1)*r
        r2 = i*r
        if i > 1
            mul!(Li,P,(disc(N,centre,r2) - disc(N,centre,r1)).*Z)
            L[1+hN*(i-1):hN*i,:] .= Li
        else
            mul!(Li,P,Z)
            L[1:hN,:] .= Li
        end
    end

    return L, rings
end

end #end module