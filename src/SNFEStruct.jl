# Contains all the structures used and define the 1D and 2D types of the problem

# Define the 2 SNFE problems types for 1D and 2D cases
abstract type Problem end
abstract type OneDim <: Problem end
abstract type TwoDim <: Problem end

include("DomainConstruct.jl")  # File containing the structure of our problem domain
include("ProbStruct.jl")       # File containing the structures of probSNFE (input and output)
include("SolveStruct.jl")      # File containing the structures of solveSNFE