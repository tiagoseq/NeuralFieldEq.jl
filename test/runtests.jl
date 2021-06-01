using Test, SNFE

#@testset "prob" begin
    
#end

@testset "solve_1D" begin
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
    prob02 = probSNFE(in02)
    in01   = Input1D(α,v,V0,L,200,T,n,I,K,S);
    prob01 = probSNFE(in01)
    in005  = Input1D(α,v,V0,L,400,T,n,I,K,S);
    prob005= probSNFE(in005)

    V02  = solveSNFE(prob02,tj)  # Solution obtained with Δx = 0.2
    V01  = solveSNFE(prob01,tj)  # Solution obtained with Δx = 0.1
    V005 = solveSNFE(prob005,tj) # Solution obtained with Δx = 0.05

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

@testset "solve_2D" begin
    # Test the convergence in time using an analytical solution obtained in
    # Pedro Lima and Evelyn Buckwar - Numerical Solution of the Neural Field Equation
    # in the two-dimensional case. Example 1, page B975
    # The exact solution is V(x,t) = exp(-t)
    using SpecialFunctions
    tol = 0.1
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
    prob002 = probSNFE(in002)
    prob001 = probSNFE(in001)

    V002 = solveSNFE(prob002,tj) # Solution obtained with Δt = 0.02
    V001 = solveSNFE(prob001,tj) # Solution obtained with Δt = 0.01

    # e_t is the absolute difference between numerical and exact solutions
    e_t = zeros(2,length(tj))
    for i = 1:length(tj)
        e_t[1,i] = abs(V002(i)[101]-exp(-tj[i]))
        e_t[2,i] = abs(V001(i)[101]-exp(-tj[i]))

        @test isapprox(2.0,e_t[1,i]./e_t[2,i],atol = tol) # e_0.02/e_0.01 ≈ 2
    end    
end