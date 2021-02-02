module SNFE

using FFTW
using Distributions: Normal, mean
using LinearAlgebra: mul!, ldiv!

# Module containing auxiliary functions
include("AuxFunctions.jl")
using .AuxFunctions

# Define the structures types for probSNFE and solveSNFE
include("SNFEStruct.jl")

export Input1D, Input2D    # Structures to wrap up the inputs to probSNFE
export probSNFE, solveSNFE # Available functions to the user

# Interface functions
include("probSNFE.jl")
include("solveSNFE1D.jl")
include("solveSNFE2D.jl")

end