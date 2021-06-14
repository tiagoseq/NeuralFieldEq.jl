using Test, NeuralFieldEq

@testset "Convergence space" begin
    # Test the convergence using an analytical bump solution
    # Laing, Troy, Gutkin, and Ermentrou - Multiple bumps in a neuronal model of working memory
    tol = 0.1
    function W(x)
        if x != 0
            return (((1-exp(-1.8*abs(x)))*(3.5/1.8) + (-1+exp(-1.52*abs(x)))*(3/1.52))*x)/abs(x) 
        else
            return 0.0
        end
    end
    a = 2.28978
    exact_solution(x) = W(x) - W(x-a)
    I(x,t) = 0.0
    K(x) = 3.5*exp(-1.8*abs(x))-3*exp(-1.52*abs(x))
    S(V) = convert(Float64,V>0.0) # Heavyside function H(V)
    α  = 1.0
    v  = 900000.0
    V0 = exact_solution
    L  = 20
    T  = 20.0
    n  = 200
    tj = [0.2,1.0,5.0,10.0,20.0]

    in02   = Input1D(α,v,V0,L,100,T,n,I,K,S);
    prob02 = probNFE(in02)
    in01   = Input1D(α,v,V0,L,200,T,n,I,K,S);
    prob01 = probNFE(in01)
    in005  = Input1D(α,v,V0,L,400,T,n,I,K,S);
    prob005= probNFE(in005)

    V02  = solveNFE(prob02,tj)  # Solution obtained with Δx = 0.2
    V01  = solveNFE(prob01,tj)  # Solution obtained with Δx = 0.1
    V005 = solveNFE(prob005,tj) # Solution obtained with Δx = 0.05

    # e_x is the absolute difference between numerical and exact solutions
    e_x = zeros(3,length(tj))
    for i = 1:length(tj)
        e_x[1,i] = abs(V02(i)[51]-exact_solution(0))
        e_x[2,i] = abs(V01(i)[101]-exact_solution(0))
        e_x[3,i] = abs(V005(i)[201]-exact_solution(0))

        @test isapprox(2.0,e_x[1,i]./e_x[2,i],atol = tol) # e_0.2/e_0.1  ≈ 2
        @test isapprox(2.0,e_x[2,i]./e_x[3,i],atol = tol) # e_0.1/e_0.05 ≈ 2
    end
end

@testset "Convergence time" begin
    # Test the convergence in time using an analytical solution obtained in
    # Pedro Lima and Evelyn Buckwar - Numerical Solution of the Neural Field Equation
    # in the two-dimensional case. Example 1, page B975
    # The exact solution is V(x,t) = exp(-t)
    tol = 0.1

    function erf(x)
        # Method to compute the error function
    
        signal = sign(x) # save the sign of x
        x = abs(x)
    
        # constants
        a1 =  0.254829592
        a2 = -0.284496736
        a3 =  1.421413741
        a4 = -1.453152027
        a5 =  1.061405429
        p  =  0.3275911
    
        t = 1.0/(1.0 + p*x)
        y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*exp(-x*x)
        return signal*y # erf(-x) = -erf(x)
    end

    function I(x,y,t)
        λ = 8.0
        b = pi/λ*(erf(sqrt(λ)))^2 #-1 to 1
        return -tanh(exp(-t))*b
    end
    function K(x,y)
        λ = 8.0
        return exp(-λ*(x^2+y^2))
    end
    S(V) = tanh(V)
    
    α  = 1.0
    v  = 900000.0
    V0 = 1.0
    L  = 2
    N  = 200
    T  = 0.1
    tj = [0.02,0.04,0.06,0.08,0.1]
    
    in002   = Input2D(α,v,V0,L,N,T,5,I,K,S);
    in001   = Input2D(α,v,V0,L,N,T,10,I,K,S);
    prob002 = probNFE(in002)
    prob001 = probNFE(in001)

    V002 = solveNFE(prob002,tj) # Solution obtained with Δt = 0.02
    V001 = solveNFE(prob001,tj) # Solution obtained with Δt = 0.01

    # e_t is the absolute difference between numerical and exact solutions
    e_t = zeros(2,length(tj))
    for i = 1:length(tj)
        e_t[1,i] = abs(V002(i)[101]-exp(-tj[i]))
        e_t[2,i] = abs(V001(i)[101]-exp(-tj[i]))

        @test isapprox(2.0,e_t[1,i]./e_t[2,i],atol = tol) # e_0.02/e_0.01 ≈ 2
    end    
end

@testset "Delayed NFE 2D" begin
    # Neural Field called breather
    # Example inspired in A. Hutt and N. Rougier - Numerical simulation scheme of one- 
    # and two-dimensional neural fields involving space-dependent delays
    function K(x,y)
        A = 20.0/(10.0*pi)
        B = 14.0/(18.0*pi)
        
        return A*exp(-sqrt(x^2+y^2)) - B*exp(-sqrt(x^2+y^2)/3.0)
    end
    S(V)=convert(Float64,V>0.005) # Heavyside function H(V-Vthresh)
    I(x,y,t)=(5.0/(32.0*pi))*exp(-(x^2+y^2)/32.0)

    α  = 1.0
    v  = 20.0
    V0 = 0.0
    L  = 20
    N  = 256
    T  = 40.0
    n  = 400
    tj = collect(0:0.2:T);
    in2d = Input2D(α,v,V0,L,N,T,n,I,K,S);

    prob = probNFE(in2d)
    V    = solveNFE(prob,tj)
    
    # Test calable structures deterministic 2D
    @test maximum(V(40.0)) == maximum(V(201))
end

@testset "Delayed SNFE 2D" begin
    # Neural Field called breather
    # Example inspired in A. Hutt and N. Rougier - Numerical simulation scheme of one- 
    # and two-dimensional neural fields involving space-dependent delays
    function K(x,y)
        A = 20.0/(10.0*pi)
        B = 14.0/(18.0*pi)
         
        return A*exp(-sqrt(x^2+y^2)) - B*exp(-sqrt(x^2+y^2)/3.0)
    end
    S(V)=convert(Float64,V>0.005) # Heavyside function H(V-Vthresh)
    I(x,y,t)=(5.0/(32.0*pi))*exp(-(x^2+y^2)/32.0)
 
    α  = 1.0
    v  = 20.0
    V0 = 0.0
    L  = 20
    N  = 256
    T  = 40.0
    n  = 400
    tj = collect(0:0.2:T);
    in2d = Input2D(α,v,V0,L,N,T,n,I,K,S);
 
    prob = probNFE(in2d)
    Vsto = solveNFE(prob,tj,0.0005,10) # 10 trajectories, ϵ=0.0005, ξ=0.1 (default value)
     
    # Test calable structures stochastic 2D
    # Mean sample solution
    @test maximum(Vsto(40.0)) == maximum(Vsto(201))

    # Trajectory 4
    @test maximum(Vsto(40.0,4)) == maximum(Vsto(201,4))
end

@testset "Delayed NFE 1D" begin
    # Example in Kulikov, Kulikova and Lima - Numerical simulation of Neural Fields
    # with Finite Transmission Speed and Random Disturbance
    function I(x,t)
        σ  = 3.0   # External input parameter
        if t <= 1.0
            return -2.89967 + 8.0*exp(-x^2/(2.0*σ^2)) - 0.5
        else
            return -2.89967
        end
    end
    K(x) = 2*exp(-0.08*sqrt(x^2))*(0.08*sin(pi*sqrt(x^2)/10)+cos(pi*sqrt(x^2)/10))
    S(V) = convert(Float64,V>0.0) # Heavyside function H(V)
    
    α  = 1.0
    v  = 25.0
    V0 = 0.0
    L  = 200
    N  = 1000
    T  = 10.0
    n  = 500
    in1D = Input1D(α,v,V0,L,N,T,n,I,K,S);

    prob = probNFE(in1D)
    V    = solveNFE(prob,[1.0,5.0,10.0])

    # Test calable structures deterministic 1D
    @test maximum(V(10.0)) == maximum(V(3))
end
 
@testset "Delayed SNFE 1D" begin
    # Example in Kulikov, Kulikova and Lima - Numerical simulation of Neural Fields
    # with Finite Transmission Speed and Random Disturbance
    function I(x,t)
        σ  = 3.0   # External input parameter
        if t <= 1.0
            return -2.89967 + 8.0*exp(-x^2/(2.0*σ^2)) - 0.5
        else
            return -2.89967
        end
    end
    K(x) = 2*exp(-0.08*sqrt(x^2))*(0.08*sin(pi*sqrt(x^2)/10)+cos(pi*sqrt(x^2)/10))
    S(V) = convert(Float64,V>0.0) # Heavyside function H(V)
    
    α  = 1.0
    v  = 15.0
    V0 = 0.0
    L  = 200
    N  = 1000
    T  = 10.0
    n  = 500
    in1D = Input1D(α,v,V0,L,N,T,n,I,K,S);

    prob = probNFE(in1D)
    Vsto = solveNFE(prob,[1.0,5.0,10.0],0.01,10) # 10 trajectories, ϵ=0.01, ξ=0.1 (default value)

    # Test calable structures stochastic 1D
    # Mean sample solution
    @test maximum(Vsto(10.0)) == maximum(Vsto(3))

    # Trajectory 4
    @test maximum(Vsto(10.0,4)) == maximum(Vsto(3,4))
end