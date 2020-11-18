module Peel

using FFTW: fft, ifft, fftshift, ifftshift
export peel

function disc(N,centre,radius)
    # This function constructs a disc in a matrix
    # 0 0 1 0 0
    # 0 1 0 1 0
    # 0 0 1 0 0

    D = Array{Float64,2}(undef,N,N)
    D = [sqrt((i-centre[1])^2+(j-centre[2])^2)<=radius ? D[i,j] = 1.0 : D[i,j] = 0.0 for i = 1:N, j = 1:N]

    return D

end

function peel(Z,centre,r,N)
    # Compute the maximum diameter to get number of rings
    rx = max(N-centre[1]-1,centre[1]-1)
    ry = max(N-centre[2]-1,centre[2]-1)
    R  = sqrt(rx^2+ry^2)

    rings = 1+floor(Int64, R/r) # Number of delay rings
    L = Array{Complex{Float64},2}(undef,N*rings,N)

    for i = 1:rings
        r1 = (i-1)*r
        r2 = i*r
        if i > 1
            L[1+N*(i-1):N*i,:] = fftshift(fft(ifftshift((disc(N,centre,r2) - disc(N,centre,r1)).*Z)))
        else
            L[1:N,:] = (disc(N,centre,r2) - disc(N,centre,r1)).*Z
        end
    end

    L[centre[1],centre[2]] = Z[centre[1],centre[2]]
    L[1:N,:] = fftshift(fft(ifftshift(L[1:N,:])))

    return L, rings
end

end #end module