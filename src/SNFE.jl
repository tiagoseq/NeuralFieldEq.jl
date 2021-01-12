module SNFE

using FFTW, Distributions, Statistics, LinearAlgebra

# Module containing auxiliary functions
include("AuxFunctions.jl")
using .AuxFunctions

export probSNFE, solveSNFE
export Input1D, Input2D # Input structures

# Define the structures types for probSNFE and solveSNFE
include("SNFEStruct.jl")
#include("DiscretiseDomain.jl")
#include("ProbStruct.jl")
#include("SolveOutput.jl")

# Interface functions
include("probSNFE.jl")
include("solveSNFE1D.jl")
include("solveSNFE2D.jl")


end