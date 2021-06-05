# Neural Field Equations

[![Build Status](https://travis-ci.com/tiagoseq/NFE.jl.svg?branch=master)](https://travis-ci.com/tiagoseq/NFE.jl)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/tiagoseq/NFE.jl?svg=true)](https://ci.appveyor.com/project/tiagoseq/NFE-jl)
[![Coverage](https://codecov.io/gh/tiagoseq/NFE.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/tiagoseq/NFE.jl)
[![Coverage](https://coveralls.io/repos/github/tiagoseq/NFE.jl/badge.svg?branch=master)](https://coveralls.io/github/tiagoseq/NFE.jl?branch=master)

Deterministic equation with delay:

<a href="https://www.codecogs.com/eqnedit.php?latex=\alpha&space;\frac{\partial&space;V}{\partial&space;t}\left(\mathbf{x},t\right)&space;=&space;I\left(\mathbf{x},t\right)&space;-&space;V\left(\mathbf{x},t\right)&space;&plus;&space;\int_{\Omega}&space;K\left(||\mathbf{x}-\mathbf{y}||_2\right)S\big[V\left(\mathbf{y},t-d\left(\mathbf{x},\mathbf{y}\right)\right)\big]\,\,d^2\mathbf{y}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\alpha&space;\frac{\partial&space;V}{\partial&space;t}\left(\mathbf{x},t\right)&space;=&space;I\left(\mathbf{x},t\right)&space;-&space;V\left(\mathbf{x},t\right)&space;&plus;&space;\int_{\Omega}&space;K\left(||\mathbf{x}-\mathbf{y}||_2\right)S\big[V\left(\mathbf{y},t-d\left(\mathbf{x},\mathbf{y}\right)\right)\big]\,\,d^2\mathbf{y}" title="\alpha \frac{\partial V}{\partial t}\left(\mathbf{x},t\right) = I\left(\mathbf{x},t\right) - V\left(\mathbf{x},t\right) + \int_{\Omega} K\left(||\mathbf{x}-\mathbf{y}||_2\right)S\big[V\left(\mathbf{y},t-d\left(\mathbf{x},\mathbf{y}\right)\right)\big]\,\,d^2\mathbf{y}" /></a>

Stochastic equation with delay:

![image](https://user-images.githubusercontent.com/73851660/120825057-eb42a100-c550-11eb-9f91-c4b4401744a9.png)
# Numerical solver for deterministic or stochastic NFE with or without delay in 1D and 2D

A numerical method was developed to solve SNFEs within the scope of the thesis "Numerical simulations of one- and two-dimensional stochasticneural field equations with delay", under the supervision of professor Pedro M. Lima.
In the deterministic case it was implemented the novel numerical method published originally by A. Hutt and N. Rougier [1,2]. Regarding the stochastic case, the approach taken was to combine the mentioned work [1,2] with the Galerkin-type method developed by C. Kuehn and M. Riedler [3], where the authors derive a spectral method to handle with stochastic neural fields with additive noise, <a href="https://www.codecogs.com/eqnedit.php?latex=\inline&space;\large&space;\epsilon" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\inline&space;\large&space;\epsilon" title="\large \epsilon" /></a>, and spatial correlation, <a href="https://www.codecogs.com/eqnedit.php?latex=\inline&space;\large&space;\xi" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\inline&space;\large&space;\xi" title="\large \xi" /></a>.

With respect to time discretisation, the scheme implemented is the explicit Euler (deterministic case) and the analogous Euler-Maruyama (stochastic case). More numerical are to be added, such as the second order scheme presented in [4] (deterministic case) and the Milstein method and the Itô-Taylor Expansion (stochastic case). In space the Fast Fourier Transforms were used to compute the solution.

[1] A. Hutt and N. Rougier - "Numerical simulation scheme of one and two-dimensional neural fields involving space-dependent delays"

[2] A. Hutt and N. Rougier - "Neural Field Simulator: two-dimensional spatio-temporal dynamics involving finite transmission speed"

[3] C. Kuehn and M.Riedler - "Large Deviations for Nonlocal Stochastic Neural Fields"

[4] P. Lima and E. Buckwar - "Numerical Solution of the Neural Field Equation in the Two-Dimensional Case"

## Dependencies

This package needs the previous installation of the packages `FFTW.jl`, `Distributions.jl` and `LinearAlgebra.jl`.

## Installation

To install the package, type on the Julia REPL the following:
```julia
using Pkg
Pkg.add("NeuralFieldEq")
```

## Solving neural field equations

Some examples are shown below in order to exemplify the correct usage of the numerical solver.

The following example shows an one-dimensional neural field:
```julia
using NeuralFieldEq # Pre-compile the package

### Defining the inputs
I(x,t) = t<=1 ? -2.89967 + 8.0*exp(-x^2/(2.0*3^2)) - 0.5 : -2.89967              # External input
K(x)   = 2*exp(-0.08*sqrt(x^2))*(0.08*sin(pi*sqrt(x^2)/10)+cos(pi*sqrt(x^2)/10)) # Connectivity function
S(V)   = convert(Float64,V>0.0) # Heavyside function H(V)                        # Firing rate. Heavyside function H(V)

α  = 1.0   # Constant decay      (must be float)
v  = 20.0  # Finite axonal speed (must be float)
V0 = 0.0   # Initial condition (can be a constant or a function)
L  = 100   # Domain length     (can be a integer or float)
N  = 512   # Number of nodes to discretise space (must be integer)
T  = 20.0  # Time span (must be float)
n  = 200   # Number of nodes to discretise time  (must be integer)
###

### Prepare the data to be solved
in_1D   = Input1D(α,v,V0,L,N,T,n,I,K,S); # Wrap the inputs in structure Input1D
prob_1D = probNFE(in_1D) # Prepare the data to solve the NFE

### Compute the deterministic solution
tj   = [5.0,10.0,20.0]      # Choose the time instants where the solution is saved
V_1D = solveNFE(prob_1D,tj) # Solve the equation and save at tj instants
```

The stochastic version of the same neural field can be computed by 
```julia
# Solve the stochastic equation with ϵ = 0.05 and ξ=0.1
# for 100 trajectories and save at tj instants
Vsto_1D = solveNFE(prob_1D,tj,0.05,100) # ξ default value is 0.1 (spatial correlation)

# Choosing another value for ξ 
Vsto2_1D = solveNFE(prob_1D,tj,0.05,100,0.15) # ξ = 0.15
```

A 2D version of the latter example can be setted using the structure `Input2D` and the functions I, K and S with the proper arguments:
```julia
# In 2D we have 2 spacial variables, x and y
I(x,y,t) = t<=1 ? -2.89967 + 8.0*exp(-x^2-y^2/(2.0*3^2)) - 0.5 : -2.89967
K(x,y)   = 2*exp(-0.08*sqrt(x^2+y^2))*(0.08*sin(pi*sqrt(x^2+y^2)/10)+cos(pi*sqrt(x^2+y^2)/10))


in_2D   = Input1D(α,v,V0,L,N,T,n,I,K,S); # Wrap the inputs in structure Input2D
prob_2D = probNFE(in_2D) # Prepare the data to solve the NFE

### Compute the two-dimensional solutions
V_2D    = solveNFE(prob_2D,tj) # deterministic
Vsto_2D = solveNFE(prob_2D,tj,0.05,100) # stochastic
```
## License
[MIT](https://choosealicense.com/licenses/mit/)
