module SNFE

import Base: show

using FFTW
using Distributions: mean
using LinearAlgebra: mul!, ldiv!

# Module containing auxiliary functions
include("AuxFunctions.jl")
using .AuxFunctions

# Define the structures types for probSNFE and solveSNFE
include("SNFEStruct.jl")

export Input1D, Input2D    # Structures to wrap up the inputs to probSNFE
export probSNFE, solveSNFE, SolveOutSto1D, SolveOutSto2D # Available functions to the user

# Interface functions
include("probSNFE.jl")
include("solveSNFE1D.jl")
include("solveSNFE2D.jl")
include("printOut.jl") # Printing outputs

end