module SNFE

using FFTW, Distributions

# Include module containing several auxiliary functions
include("Peel.jl")
using .Peel

export probSNFE, solveSNFE

# Include data structures
include("DiscretiseDomain.jl")
include("ProbOutput.jl")
include("SolOutput.jl")

# Include solver functions
include("probSNFE.jl")
include("solveSNFE.jl")


end