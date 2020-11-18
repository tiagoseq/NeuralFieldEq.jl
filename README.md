# SNFE

[![Build Status](https://travis-ci.com/tiagoseq/SNFE.jl.svg?branch=master)](https://travis-ci.com/tiagoseq/SNFE.jl)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/tiagoseq/SNFE.jl?svg=true)](https://ci.appveyor.com/project/tiagoseq/SNFE-jl)
[![Coverage](https://codecov.io/gh/tiagoseq/SNFE.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/tiagoseq/SNFE.jl)
[![Coverage](https://coveralls.io/repos/github/tiagoseq/SNFE.jl/badge.svg?branch=master)](https://coveralls.io/github/tiagoseq/SNFE.jl?branch=master)


# Stochastic Neural-Field Equations

A numerical method written in julia to solve Stochastic Neural-Field Equations within the scope of my master's thesis. This method was published by Axel Hutt and Nicolas P. Rougier in "Numerical simulation scheme of one and two-dimensional neural fields involving space-dependent delays".


## Installation

To install and pre-compile the package type:

```julia
Pkg.add("SNFE")
using SNFE
```

## Dependencies

This package needs the previous installation of the packages FFTW.jl and Distributions.jl

## Usage

```julia
using SNFE
prob = probSNFE(inputs)
sol  = solveSNFE(prob)
```

## License
[MIT](https://choosealicense.com/licenses/mit/)