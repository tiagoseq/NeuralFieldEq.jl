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

export Input1D, Input2D    # Public structures to wrap up the inputs to probSNFE
export probSNFE, solveSNFE # Public functions to compute SNFE
export SolveOutSto1D, SolveOutSto2D # Structures to print output of solveSNFE

# Interface functions
include("probSNFE.jl")
include("solveSNFE1D.jl")
include("solveSNFE2D.jl")

# Format outputs for probSNFE and solveSNFE
include("printOut.jl")

end