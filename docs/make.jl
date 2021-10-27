push!(LOAD_PATH,"../src/")

using NeuralFieldEq, Documenter

makedocs(
    sitename = "NeuralFieldEq.jl",
    modules  = [NeuralFieldEq],
    authors  = "Tiago Sequeira",
    pages    = [
                "Home" => "index.md"
                "User guide" => "usage.md"
                "Examples" => "examples.md"
                "API documentation" => "apidocs.md"
               ])

deploydocs(repo="github.com/tiagoseq/NeuralFieldEq.jl",)