module NeuralFieldEq

import Base: show

using FFTW
using Distributions: mean
using LinearAlgebra: mul!, ldiv!
using ProgressMeter: Progress, next!

# Module containing auxiliary functions
include("AuxFunctions.jl")
using .AuxFunctions

# Define the structures types for probNFE and solveNFE
include("NFEStruct.jl")

export Input1D, Input2D  # Public structures to wrap up the inputs to probNFE
export probNFE, solveNFE # Public functions to compute NFE
export SolveOutSto1D, SolveOutSto2D # Structures to print output of solveNFE

# Interface functions
include("probNFE.jl")
include("solveNFE1D.jl")
include("solveNFE2D.jl")

# Format outputs for probSNFE and solveSNFE
include("printOut.jl")

end