function show(io::IO, p::ProbOutput1D)
    print(io,
        " ├─ Domain:       Ω × [0,T] = [$(p.Ω.x[1]),$(p.Ω.x[end])] × [0,$(p.Ω.t[end])]", '\n',
        " ├─ Spatial step: dx   = $(p.Ω.dx)", '\n',
        " ├─ Time step:    dt   = $(p.Ω.dt)", '\n',
        " ├─ Velocity:     v    = $(p.Ω.v)", '\n',
        " ├─ Delay rings:  umax = $(p.Ω.rings1D)")
end

function show(io::IO, p::ProbOutput2D)
    print(io,
        " ├─ Domain:       Ω × [0,T] = [$(p.Ω.x[1]),$(p.Ω.x[end])]^2 × [0,$(p.Ω.t[end])]", '\n',
        " ├─ Spatial step: dx   = $(p.Ω.dx)", '\n',
        " ├─ Time step:    dt   = $(p.Ω.dt)", '\n',
        " ├─ Velocity:     v    = $(p.Ω.v)", '\n',
        " ├─ Delay rings:  umax = $(p.Ω.rings2D)")
end

function show(io::IO, s::SolveOutDet1D)
    if issubset(s.tsaved,s.t)
        print(io,
              " ├─ Solution saved at all time instants.", '\n',
              " ├─ 1D deterministic solution")
    else
        print(io,"Solution not saved at all time instants. Choose appropriate time instants.")
    end
end
function show(io::IO, s::SolveOutDet2D)
    if issubset(s.tsaved,s.t)
        print(io,
        " ├─ Solution saved at all time instants.", '\n',
        " ├─ 2D deterministic solution")
    else
        print(io,"Solution not saved at all time instants. Choose appropriate time instants.")
    end
end
function show(io::IO, s::SolveOutSto1D)
    if issubset(s.tsaved,s.t)
        print(io,
        " ├─ Solution saved at all time instants.", '\n',
        " ├─ 1D stochastic solution")
    else
        print(io,"Solution not saved at all time instants. Choose appropriate time instants.")
    end
end
function show(io::IO, s::SolveOutSto2D)
    if issubset(s.tsaved,s.t)
        print(io,
        " ├─ Solution saved at all time instants.", '\n',
        " ├─ 2D stochastic solution")
    else
        print(io,"Solution not saved at all time instants. Choose appropriate time instants.")
    end
end