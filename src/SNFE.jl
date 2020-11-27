module SNFE

using FFTW, Distributions, Statistics, LinearAlgebra

# Module containing several auxiliary functions
include("AuxFunctions.jl")
using .AuxFunctions

export probSNFE, solveSNFE

# Include data structures
include("DiscretiseDomain.jl")
include("ProbOutput.jl")
include("SolOutput.jl")

# Include solver functions
include("probSNFE.jl")
include("solveSNFE.jl")


end