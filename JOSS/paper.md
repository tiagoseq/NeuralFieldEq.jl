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
 - name: Instituto Superior Técnico - University of Lisbon
   index: 1
date: 02 November 2021
bibliography: paper.bib
---

# Summary

In their works, @WilsonCowan presented a model, later developed by @Amari, named Neural Field Equations (NFE). @Bressloff:minibook and @Coombes:WavesBumps, summarise the main theoretical results and enumerate a wide variety of neurological phenomena that can be described by these equations. Such as working memory [@Laing:WM], travelling waves, etc. One of the improvements to consider in the simpler NFE, is to take into account the time taken by the stimulus to travel from neurons at different points. This leads to the delayed equation,
\begin{equation}\label{eq:dNFE}
  \alpha \frac{\partial V}{\partial t}\left(\mathbf{x},t\right) = I\left(\mathbf{x},t\right) - V\left(\mathbf{x},t\right) + \int_{\Omega} K\left(||\mathbf{x}-\mathbf{y}||_2\right)S\big[V\left(\mathbf{y},t-d\left(\mathbf{x},\mathbf{y}\right)\right)\big]\,\,d^k\mathbf{y},
\end{equation}
with $\Omega=[-\frac{L}{2},\frac{L}{2}]^k$ with $k=1,2$; $V(\mathbf{x},t)$ is the membrane potential at point $\mathbf{x} \in \Omega$ at time $t$; $I(\mathbf{x},t)$ is the external input applied to the neural field; $K\left(||\mathbf{x}-\mathbf{y}||_2\right)$ is the average strength of connectivity between neurons located at points $\mathbf{x}$ and $\mathbf{y}$ in an isotropic and homogeneous neural field. When the coupling is positive (resp. negative) the synapses are excitatory (resp. inhibitory); $S(V)$ is the firing rate function, which converts the potential to the respective firing rate result; $\alpha$ is the constant decay rate; $d\left(\mathbf{x},\mathbf{y}\right)=\frac{||\mathbf{x}-\mathbf{y}||_2}{v}$ represents the delay, which is assumed to only depend on the distance, between points $\mathbf{x}$ and $\mathbf{y}$, and on the transmission speed $v$. If $v$ is sufficiently high, $d$ can be neglected and \autoref{eq:dNFE} is reduced to a non-delayed NFE. Moreover, $V$ satisfies the initial condition $V(\mathbf{x},t_0)=V_0(\mathbf{x},t_0),\,\,t_0\in\left[-\tau_{max},\,\,0\right]$, wherein $\tau_{max}$ is the maximum delay correspondent to $\Omega$.

@KuehnRiedler considered the stochasticity inherent to neuronal interactions, and presented a spectral method to stochastic non-delayed NFE with additive white noise and spatial correlation. So, considering a finite speed in this scenario, the following delayed stochastic NFE can be written,
\begin{align}\label{eq:dSNFE}
    \alpha\, dV\left(\mathbf{x},t\right) =& \left[I\left(\mathbf{x},t\right) - V\left(\mathbf{x},t\right) + \int_{\Omega}K\left(||\mathbf{x}-\mathbf{y}||_2\right)S\big[V\left(\mathbf{y},t-d\left(\mathbf{x},\mathbf{y}\right)\right)\big]\,\,d^2\mathbf{y}\right]dt + \nonumber \\
    & \epsilon dW\left(\mathbf{x},t\right),
\end{align}
whereas $\epsilon$ is the additive noise level and $W$ is a $Q$-Wiener process. With $V(\mathbf{x},t_0)=V_0(\mathbf{x},t_0),\,\,t_0\in\left[-\tau_{max},\,\,0\right]$ and $V_0$ is a given stochastic process.

# Statement of need

Studies suggest that Neural Field Equations, can be applied to cognitive robotics. The architecture of autonomous robots, able to interact with other agents in solving a mutual task, is strongly inspired by the processing principles and the neuronal circuitry in the primate brain [@ErlhagenBicho]. Furthermore, recent efforts made by @FerreiraEtAl:RapidLearning show a necessity of efficient numerical methods capable of computing NFE in real time, in order to endow robots with working memory mechanisms.

The common numerical methods for Neural Field Equations use the classical quadrature methods, which have $\mathcal{O}\left(N^{2k}\right)$ complexity, where $k$ is the dimension of the domain, therefore, they do not scale well. In the case of non-delayed equations, one can directly apply Fast Fourier Transforms (FFT), $\mathcal{O}\left(N^k\log^k N\right)$, to compute the numerical solution. However, when $v<\infty$, this is no longer the case. Motivated by this, @HuttRougier:2 proposed a novel numerical scheme, where the authors rewrote the integral in a suitable way, such that, the delayed NFE could be numerically approximated using FFT techniques.

`NeuralFieldEq.jl` aims to enhance the method derived in [@HuttRougier:2] and implemented by @Simulator, and adapt to the stochastic scenario presented by @KuehnRiedler. The solver uses Real Fast Fourier Transforms (RFFT) provided by @fftw, reducing the overall computations. And is written in `Julia` [@Julia], which can be performant like low-level languages without sacrificing the usual features of high-level languages. Furthermore, its built to handle NFE in the scenarios aforementioned and their respective combinations, i.e. non-delayed or delayed equations in deterministic or stochastic 1D or 2D neural fields.

@SequeiraLima exploited this versatility and efficiency to study stochastic neural fields where low transmission velocities play an important role, facilitating the exploration of new ideas and experiments.

# Package usage

1. Specify parameters, initial condition and functions by using `Input1D` or `Input2D`:
    - The functions are defined as: External input: `I(x,t)` or `I(x,y,t)`; Kernel: `K(x)` or `K(x,y)`; And Firing rate: `S(V)`.
    - The initial condition can be constant or a function. In the latter case, must be defined as `V0(x)` or `V0(x,y)`;
    - The parameters are: `a`-$\alpha$, `v`-velocity, `L`-domain length, `N`-number of spatial nodes, `T`-time span and `n`-number of time steps;
    - The inputs must be wrapped in the following order: `nf = Input1D(a,v,V0,L,N,T,n,I,K,S)`;
    - **Remark:** Currently, to obtain the non-delayed problem, the velocity must satisfy: $v>\frac{L}{\sqrt{2}\Delta t}$ in 2D and $v>\frac{L}{2\Delta t}$ in 1D.

2. Prepare the NFE through function `probNFE`. Example: `NFE = probNFE(nf)`;

3. Solve the equation using `solveNFE`:
    - Declare an array with the time instants where the solution is saved;
    - Compute deterministic solution: `Vdet = solveNFE(NFE,tj)`;
    - Compute stochastic solution with noise magnitude `eps`, spatial correlation `xi`, for `np` paths: `Vsto = solveNFE(NFE,tj,eps,xi,np)`;
    - **Remark:** Currently, `xi=0.1` is the default value.

4. Handle the solution:
    - The output of `solveNFE` is a callable object:
    - Deterministic solution at a time instant. Example: `Vdet(t1)`;
    - $p^{th}$ trajectory at a time instant. Example: `Vsto(t1,p)`;
    - Mean stochastic solution at a time instant. Example: `Vsto(t1)`;
    - Return spatial vector: `Vdet.x`.

# Example and code performance

The example shown below is adapted from @Kulikov1D
```julia
using NeuralFieldEq, Plots
I(x,t) = -3.4 + 8.0*exp(-x^2/(2.0*3^2))
K(x) = 2*exp(-0.08*sqrt(x^2))*(0.08*sin(pi*sqrt(x^2)/10)+cos(pi*sqrt(x^2)/10))
S(V) = V<=0.0 ? 0.0 : 1.0 # Heaviside function
nf  = Input1D(1.0,20.0,0.0,100,512,20.0,200,I,K,S);
NFE = probNFE(nf)
ti  = [5.0,10.0,20.0]  # Choose instants where the sol is saved
V   = solveNFE(NFE,ti) # Compute deterministic solution
Vsto = solveNFE(NFE,ti,0.05,100) # eps = 0.05. xi = 0.1. 100 paths
plot(V.x,[V(1),Vsto(1),Vsto(1,4)],xlabel="x",ylabel="Action Potential",
     label=["Det Solution" "Sto Mean Solution" "4th path"])
```
![Caption for example figure.\label{fig:example}](plots1D.png){width=75%}

Regarding the solver performance, the table below shows the time spent (in seconds) in computing one time step with `N` nodes, for the example listed above and its `2D` version. Both equations were computed with `v=400`.
\begin{table}[H]
\begin{tabular}{c|c|c}
$N$  & 1D     & 2D     \\ \hline
128  & 8.6e-6 & 1.4e-3 \\
256  & 2.1e-5 & 9.2e-3 \\
512  & 3.1e-5 & 3.8e-2 \\
1024 & 6.2e-5 & 0.155  \\
2048 & 1.3e-4 & 0.621  \\
4096 & 2.6e-4 & 2.72
\end{tabular}
\end{table}

# Acknowledgements

I want to thank my advisor, Pedro Lima, that kindly reviewed this article.

# References
