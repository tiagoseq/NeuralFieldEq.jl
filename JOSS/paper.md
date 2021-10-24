---
title: 'NeuralFieldEq.jl: A flexible solver to compute Neural Field Equations in several scenarios'
tags:
  - Julia
  - Computational Neuroscience
  - Neural Field Equations
  - Stochastic Processes
  - Delayed Equations
  - Fast Fourier Transforms
authors:
  - name: Tiago Sequeira^[first author]
    orcid: 0000-0001-8579-3676
    affiliation: "1"
affiliations:
 - name: Instituto Superior T\'{e}cnico
   index: 1
date: 24 October 2021
bibliography: paper.bib
---

# Summary

Since @hodgkinhuxley that mathematical modelling proved to be a prominent neuroscience research area. One of its main challenges is to describe phenomena like memory, perception, etc. These neurobiological processes occur in the cortex, a layer of the brain that contains approximately $10^{10}$ neurons distributed along $2500\,cm^2$ of surface area with $2.8\,mm$ of depth. Driven by such complex structure, @WilsonCowan and @Amari derived a model, Neural Field Equation (NFE), that treats the cerebral cortex as a continuum space being able to deal with these large scale dynamical neuronal patterns.

\begin{equation}\label{eq:NFE}
  \alpha \frac{\partial V}{\partial t}\left(\mathbf{x},t\right) = I\left(\mathbf{x},t\right) - V\left(\mathbf{x},t\right) + \int K\left(||\mathbf{x}-\mathbf{y}||_2\right)S\big[V\left(\mathbf{y},t\right)\big]\,\,d^2\mathbf{y}.
\end{equation}

Since then there have been endeavours to improve NFEs, from taking into account the speed of signal propagation to considering the stochastic behaviour of the neuron, these are some improvements that aim a better description of reality.

# Statement of need

The classical quadrature methods to numerically approximate the integral present in \autoref{eq:NFE} have a complexity of $\mathcal{O}^{2}$ or $\mathcal{O}^{4}$ in 1D or 2D domains, respectively, making these methods unsuitable for efficient numerical approximations of NFEs. Although the this integral can be seen as a convolution, when considering finite signal velocities is no longer the case. @HuttRougier:2 proposed a novel numerical scheme that addresses the delayed version of \autoref{eq:NFE} that rewrites the integral into a convolutional form in order to apply a Fourier Transform to it, implying a substantially speed-up when computing delayed NFEs.

The numerical method that `NeuralFieldEq.jl` implements arose from the combination of the key idea developed by Hutt and Rougier for delayed neural fields in the stochastic scenario presented by @RiedlerChristian, where the authors proved the convergence of spectral methods applied to stochastic NFEs with additive white noise spatially correlated.

`Julia` [@Julia] code when well written, designed and profiled its performance can be close to `C` or `Fortran` without sacrificing the usual features present in high-level languages. Also, the package makes use of the the multiple dispatch concept, allowing it to be flexible enough to handle NFEs in three different scenarios, 1D or 2D domains, non-delayed or delayed equations and deterministic or stochastic neural fields. These advantages were the trigger needed to develop a new user friendly and fast NFE solver, improving the `Python` solver written by @Simulator]

# Package usage

To illustrate the code usage we will discuss an one-dimensional example taken from @Kulikov1D.

The first step is the initialisation of inputs followed by wrapping it in structures `Input1D` or `Input2D`, depending on the domain dimension. A remark: The input function has to have `x`,`t` or `x`,`y`,`t` as its arguments, the connectivity function with `x` or `x`,`y` and the firing rate with `V`.
```julia
using NeuralFieldEq
# 1D neural field
# Define function inputs, I, K and S
I(x,t) = -2.89967 + 8.0*exp(-x^2/(2.0*3^2)) - 0.5
K(x) = 2*exp(-0.08*sqrt(x^2))*(0.08*sin(pi*sqrt(x^2)/10)+cos(pi*sqrt(x^2)/10))
S(V) = V<=1.0 ? 0.0 : 1.0

# Define parameters
a  = 1.0  # Constant decay      
v  = 20.0 # Finite axonal speed
V0 = 0.0  # Initial condition
L  = 100  # Domain length
N  = 512  # Number of nodes to discretise space
T  = 20.0 # Time span
n  = 200  # Number of nodes to discretise time

nf_1d = Input1D(a,v,V0,L,N,T,n,I,K,S); # Wrap inputs in structure Input1D
```
The second step is to use the function `probNFE` to prepare the inputs to solve the equation and finally computing the solution using `solveNFE`.
```julia
tj   = [5.0,10.0,20.0]   # Choose instants where the sol is saved
prob = probNFE(nf_1d)    # Prepare the data to solve the NFE
V    = solveNFE(prob,tj) # Solve the equation and save at tj
```
If we want to address the stochastic case we simply need to type extra arguments in `solveNFE`.
```julia
# Solve the stochastic equation 100 times
# Noise magnitude: eps = 0.05. Correlation coefficient: xi = 0.1
Vsto  = solveNFE(prob,tj,0.05,100)      # xi default value 0.1
Vsto2 = solveNFE(prob,tj,0.05,100,0.15) # xi = 0.15
```
The variables obtained with `solveNFE` are equipped with practical features to handle the solutions and create plots.
```julia
using Plots
V(10.0)     # Returns the deterministic solution at t=10.0
Vsto(20.0)  # Returns the mean stochastic solution at t=20.0
Vsto(5.0,4) # Returns the 4th trajectory at t=5.0

x = V.x # Returns the spatial vector
plot(x,[V_1D(1),Vsto_1D(1),Vsto_1D(1,4)],
     xlabel="x",
     ylabel="Action potential",
     label=["Deterministic solution"
            "Stochastic mean solution"
            "4th trajectory"])
```
![Caption for example figure.\label{fig:example}](plots.png){width=80%}

# Acknowledgements

I want to thanks my professor Pedro Lima that kindly reviewed this article.

# References
