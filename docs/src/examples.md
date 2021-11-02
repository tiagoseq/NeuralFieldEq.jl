# Examples

In this section we will address some examples to show the versatility of the package. From setting the neural field to preparing the equation to be computed, as well handling the solutions in order to facilitate the study and interpretation of it.

## How to use it
The solver is divided into three steps:
- Introduce the parameters and functions using the structures `Input1D` or `Input2D`;
- Pre-process the NFE using the function `probNFE`;
- Solve the equation using the function `solveNFE` at time instants chosen by the user, with or without noise.

With respect to the structures `Input1D` and `Input2D`, they are needed to wrap the inputs needed to define our NFE. They have the following order:
- `α  :: AbstractFloat`: Decay rate
- `v  :: AbstractFloat`: Axonal speed
- `V0 :: fV0`          : Initial condition (constant or a function)
- `L  :: Number`       : Domain's length
- `N  :: Integer`      : Number of spatial nodes
- `T  :: AbstractFloat`: Time span
- `n  :: Integer`      : Number of time nodes
- `I  :: fI`           : External input function
- `K  :: fK`           : Connectivity function
- `S  :: fS`           : Firing rate function

**Remark 1**: The function `I`, depending on the dimensionality of the domain, has to have `x`,`t` or `x`,`y`,`t` as its arguments. Function `K` has `x` or `x`,`y`. Function `S` with `V`.

Once we define our input structure, we can now pass as input to function `probNFE`, where the NFE is prepared to be solved using the function `solveNFE`.

**Remark 2**: Currently, to work with the non-delayed problem, the velocity to insert must satisfy the condition: ``v>\frac{L}{\sqrt{2}\Delta t}`` in 2D and ``v>\frac{L}{2\Delta t}`` in 1D, meaning that in practice the user must specify a big velocity (ex.: ``999999.0``).

Below is shown an abstract example of the `NeuralFieldEq.jl` usage:
```julia
nfe  = Input1D(α,v,V0,L,N,T,n,I,K,S); # Wrap the inputs in structure Input1D
prob = probNFE(nfe)                   # Pre-process the NFE to be computed

# Solve the deterministic 1D problem
Vdet = solveNFE(prob,[ti,tj,tk]) # solution saved at ti, tj, and tk

# Solve the stochastic 1D problem np times
# ϵ level of additive noise, spatial correlation (0.1 default value)
Vsto = solveNFE(prob,[ti,tj,tk],ϵ,np,xi=0.1) # solution saved at ti, tj, and tk
```

## Handling solutions

Considering `t=[ti,tj,tk]`, to access the deterministic solution at `tj`, the user must type `V(tj)` or `V(2)` (the index of `tj` is 2). In the stochastic case, the procedure is the same. Note that if the user only specifies the time instant, for example `Vsto(tk) # == Vsto(3)`, then it will be returned the **mean solution** at `tk`. Whereas to choose a **trajectory** `p` at time `tk`, the command should be `Vsto(tk,p) # == Vsto(3,p)` for the `p` trajectory at `tk`.
In the 1D case the output returned is a vector, while in a 2D space is a matrix.

The user can also obtain a specific point in space and time. So, considering the discretised space vector `x=[x1,x2,...,xN]`, the solution at point `x2` at time `tj` is accessed by typing `V(x2,tj)` (and in 2D: `V(x2,y1,tj)`). Regarding dealing with the stochastic solution, the same idea explained above still holds, `Vsto(x2,tj)` stands for the value of the **mean solution** at `(x2,tj)`, and `Vsto(x2,tj,p)` stands for the value of trajectory `p` at point `(x2,tj)`.

Moreover, to help plotting the solutions, the output of the function `solveNFE` is endowed with the fields `x`, `y` (in 2D case), and `tsaved` (instants where the solution was saved), corresponding to the discretised spatial and time variables, respectively.

## Example of deterministic and stochastic NF in 1D

This example was adapted from [Kulikov et al. (2019)]
```julia
using NeuralFieldEq
# 1D neural field
# Define function inputs, I, K and S
I(x,t) = -2.89967 + 8.0*exp(-x^2/(2.0*3^2)) - 0.5
K(x) = 2*exp(-0.08*sqrt(x^2))*(0.08*sin(pi*sqrt(x^2)/10)+cos(pi*sqrt(x^2)/10))
S(V) = V<=0.0 ? 0.0 : 1.0 # Heaviside function

# Define parameters
a  = 1.0  # Constant decay
v  = 5.0  # Finite axonal speed
V0 = 0.0  # Initial condition
L  = 100  # Domain length
N  = 512  # Number of nodes to discretise space
T  = 10.0 # Time span
n  = 100  # Number of nodes to discretise time

nfe_1d = Input1D(a,v,V0,L,N,T,n,I,K,S); # Wrap inputs in structure Input1D
```
The second step is to use the function `probNFE` to prepare the inputs to solve the equation and finally computing the solution using `solveNFE`.
```julia
tj   = 0:0.2:T;          # Choose instants where the sol is saved
prob = probNFE(nfe_1d)   # Prepare the data to solve the NFE
V    = solveNFE(prob,tj) # Solve the equation and save at tj
```
Also, the function `probNFE` prints useful information such as spatial and time step as well as the number of delay rings (if ``v>\frac{L}{2\Delta t}`` the number of rings is ``1``).

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
V(1.0)      # Returns the deterministic solution at t=1.0
Vsto(10.0)  # Returns the mean stochastic solution at t=10.0
Vsto(5.0,4) # Returns the 4th trajectory at t=5.0

x = V.x # Returns the spatial vector
plot(x,[V_1D(5.0),Vsto_1D(5.0),Vsto_1D(5.0,4)],
     xlabel="x",
     ylabel="Action potential",
     label=["Deterministic solution"
            "Stochastic mean solution"
            "4th trajectory"])
```

```@raw html
<img src='../figs/plots1D.png' width=500><br>
```

## NF with finite propagation speed

Considering the previous example, we will now compare the delayed solution with the non-delayed one.
```julia
# Setting the non-delayed neural field
v_inf = 9999999.0

nfe_1d_vinf = Input1D(a,v_inf,V0,L,N,T,n,I,K,S);
prob_vinf   = probNFE(nfe_1d_vinf)
V_vinf      = solveNFE(prob_vinf,tj)

# Animate the delayed and non-delayed solutions
anim = @animate for i = 1:length(tj)
  p1 = plot(x,V(i),title="v=5. t=$(tj[i])",ylims=(-9.3,16.7))
  p2 = plot(x,V_vinf(i),title="v=inf. t=$(tj[i])",ylims=(-9.3,16.7))
  plot(p1,p2,size=(900,400))
end
gif(anim, "bumps.gif",fps=10)
```

```@raw html
<img src='../figs/bumps.gif' width=500><br>
```

## Example of deterministic and stochastic NF in 2D

This is a two-dimensional extension of the previous example, adapted from [Kulikova et al. (2019)]
```julia
# Example with the same parameters as in 1D
# In 2D we have 2 spacial variables, x and y
I_2d(x,y,t) = -2.89967 + 8.0*exp(-x^2-y^2/(2.0*3^2)) - 0.5
K_2d(x,y)   = 2*exp(-0.08*sqrt(x^2+y^2))*(0.08*sin(pi*sqrt(x^2+y^2)/10)+cos(pi*sqrt(x^2+y^2)/10))

nfe_2d  = Input2D(a,v,V0,L,N,T,n,I_2d,K_2d,S);
prob_2d = probNFE(nfe_2d)

# Compute the two-dimensional solutions
V_2d    = solveNFE(prob_2d,tj)         # deterministic
Vsto_2d = solveNFE(prob_2d,tj,0.1,20)  # stochastic

# Plotting
x = V_2d.x; # Returns the spatial vector in x direction
y = V_2d.y; # Returns the spatial vector in y direction

p1 = plot(x,y,V_2d(10.0),
          title  = "Deterministic equation at t=10",
          xlabel = "x",
          ylabel = "y",
          zlabel = "Action potential")

p2 = plot(x,y,Vsto_2d(10.0,17),
          title  = "17th path with eps=0.1 and xi=0.1 at t=10",
          xlabel = "x",
          ylabel = "y",
          zlabel = "Action potential")
plot(p1,p2,size=(1000,500))
```

```@raw html
<img src='../figs/bumps2d.png' width=500><br>
```

## Examples with V0 as a function

### One-dimensional NF

This example was taken from section 2 of [Laing et al. (2002)](http://www.math.pitt.edu/~bard/pubs/multibump.pdf) where an analytical stationary solution was derived.
```julia
using NeuralFieldEq, Plots

function W(x)
    if x != 0
        return (((1-exp(-1.8*abs(x)))*(3.5/1.8) + (-1+exp(-1.52*abs(x)))*(3/1.52))*x)/abs(x)
    else
        return 0.0
    end
end
a2   = 2.28978
sol_exact(x) = W(x) - W(x-a2)

extInput(x,t) = 0.0
kernel(x)     = 3.5*exp(-1.8*abs(x))-3*exp(-1.52*abs(x))
firingRate(V) = V<=0.0 ? 0.0 : 1.0

a  = 1.0
v  = 90000.0   # This velocity satisfies the condition v>L/(2*Δt), thus we are considering the non-delayed problem
V0 = sol_exact # The initial condition, V0, is a function
L  = 20
N  = 512
T  = 20.0
n  = 200

input = Input1D(a,v,V0,L,N,T,n,extInput,kernel,firingRate);
prob  = probNFE(input)
V     = solveNFE(prob,[0.2,1.0,5.0,10.0,20.0])         # Compute the deterministic solution

# Plotting
x = V.x
plot(x,[V(10.0),sol_exact.(x)],
     labels=["Numerical solution" "Exact solution"],
     legend=:bottomright,
     title="Stationary problem at t=10",
     xlabel="x",
     ylabel="V(x)")
```

```@raw html
<img src='../figs/stasolution1d.png' width=500><br>
```

### Two-dimensional NF
In this example we will adapt a 2D neural field taken from section 8 of [Laing et al. (2002)](http://www.math.pitt.edu/~bard/pubs/multibump.pdf), where the initial condition is not a constant, but a function of ``x`` and ``y``, more specifically in the following form:
```math
\begin{equation}
    \begin{cases}
        V_0(x,y) = 5\,\, \text{if} -4<x<5.6 \land -12<y<4 \\
        V_0(x,y) = 0\,\, \text{otherwise.}
    \end{cases}
\end{equation}
```
```julia
using NeuralFieldEq, Plots

# Define functions
I(x,y,t) = 0
W(x,y)   = exp(-0.45*sqrt(x^2+y^2))*(0.45*sin(sqrt(x^2+y^2))+cos(sqrt(x^2+y^2))) # Connectivity function
f(V)     = 2.0*exp(-0.1/((V-1.5)^2))*convert(Float64,V>1.5) # Firing rate function (an exponential times a heaviside)

# Define parameters
a  = 1.0     # Temporal decay of the synapse
v  = 99999.0 # Finite axonal speed
L  = 40      # Domain length. [-20,20]
N  = 50      # Number of nodes to discretise space
T  = 4.0     # Time span
n  = 200     # Number of nodes to discretise time
tj = [1.0,2.0,3.0,4.0]

# Initial condition
function V_0(x,y)
    if (x > -4.0 && x < 5.6) && (y > -12.0 && y < 4.0)
        return 5.0
    else
        return 0.0
    end
end

nfe  = Input2D(a,v,V_0,L,N,T,n,I,W,f);
prob = probNFE(nfe)
V    = solveNFE(prob,tj)

# Plotting
x = V.x
y = V.y
plot(x,y,V(),st=:surface,
     xlabel = "x",
     ylabel = "y",
     zlabel = "V(x,y)",
     title  = "Example of section 8 of Laing et al. (2002)")
```

```@raw html
<img src='../figs/bumps2d_laing.png' width=500><br>
```

## Breather type instabilities

The authors in [Hutt & Rougier 2013](https://hal.inria.fr/hal-00872132/document) successfully induced in a NF breather type instabilities using only the axonal speed. Here we use an adapted NFE in order to induce the same type of instabilities.

First we address a stable field (with ``v=\infty``):
```julia
using NeuralFieldEq, Plots

I(x,y,t) = (5.0/(32.0*pi))*exp(-(x^2+y^2)/32.0)
function K(x,y)
    A = 20.0/(10.0*pi)
    B = 14.0/(18.0*pi)
    return A*exp(-sqrt(x^2+y^2)) - B*exp(-sqrt(x^2+y^2)/3.0)
end
S(V) = V<=0.005 ? 0.0 : 1.0 # H(V-th)

a  = 1.0
v  = 999999.0
V0 = 0.0
L  = 20
N  = 256
T  = 10.0
n  = 200

nfe  = Input2D(a,v,V0,L,N,T,n,I,K,S);
prob = probNFE(nfe)

# Solve the neural field equation
tj = 0:0.2:T;
V  = solveNFE(prob,tj)

# Plots
minmax_V = zeros(2,length(tj))
for i = 1:length(tj)
    minmax_V[1,i] = minimum(V(i))
    minmax_V[2,i] = maximum(V(i))
end

# Animate solution
x = V.x
y = V.y
anim = @animate for i = 1:length(tj)
  p1 = plot(x,y,V(i),st=:surface,title="t=$(tj[i])",color=:balance,zlims=(-0.43,0.54))
  p2 = plot(tj,[minmax_V[1,i],minmax_V[2,i]],label=[minimum maximum],ylims=(-0.43,0.54))
  plot(p1,p2)
end
gif(anim, "vinf_breather.gif",fps=10)
```

```@raw html
<img src='../figs/vinf_breather.gif' width=500><br>
```

Next we use the finite velocities ``v=3`` and ``v=5``:
```julia
nfe_v5  = Input2D(a,5.0,V0,L,N,80.0,1600,I,K,S);
nfe_v3  = Input2D(a,3.0,V0,L,N,80.0,1600,I,K,S);
prob_v5 = probNFE(nfe_v5)
prob_v3 = probNFE(nfe_v3)

# Solve the neural field equation
V_v5 = solveNFE(prob_v5,tj)
V_v3 = solveNFE(prob_v3,tj)

# Plotting
# Animate the delayed solutions
anim = @animate for i = 1:length(tj)
  p1 = plot(x,y,V_v5(i),st=:surface,title="v=5. t=$(tj[i])",color=:balance,zlims=(-2.8,0.95))
  p2 = plot(x,y,V_v3(i),st=:surface,title="v=3. t=$(tj[i])",color=:balance,zlims=(-2.8,0.95))
  plot(p1,p2)
end
gif(anim, "v5v3_breather.gif",fps=10)
```

```@raw html
<img src='../figs/v5v3_breather.gif' width=500><br>
```

As we can see, the field oscillates in the first seconds of the simulation and lower the velocity the more time the neural field oscillates. Indeed, if we consider ``v=2`` the field is unstable for all time span considered
```julia
nfe_v2  = Input2D(a,2.0,V0,L,N,100.0,2000,I,K,S);
prob_v2 = probNFE(nfe_v2)
V_v2    = solveNFE(prob_v2,tj)

# Plotting
minmax_Vv2 = zeros(2,length(tj))
for i = 1:length(tj)
    minmax_Vv2[1,i] = minimum(V_v2(i))
    minmax_Vv2[2,i] = maximum(V_v2(i))
end

# Animate solution
anim = @animate for i = 1:length(tj)
  p1 = plot(x,y,V_v2(i),st=:surface,title="t=$(tj[i])",color=:balance,zlims=(-3.6,0.95))
  p2 = plot(tj,[minmax_Vv2[1,i],minmax_Vv2[2,i]],label=[minimum maximum],ylims=(-3.6,0.95))
  plot(p1,p2,size=(900,400))
end
gif(anim, "v2_breather.gif",fps=10)
```

```@raw html
<img src='../figs/v2_breather.gif' width=500><br>
```