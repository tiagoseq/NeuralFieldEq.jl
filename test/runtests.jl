using Test, SNFE

@testset "prob" begin
    
end

@testset "solve" begin
    tol = 1e-3

    # Test the convergence using an analytical bump solution
    # Laing, Troy, Gutkin, and Ermentrou - Multiple bumps in a neuronal model of working memory
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
    V0 = sole
    L  = 20
    T  = 20.0
    n  = 200

    input = Input1D(α,v,V0,L,N,T,n,I,K,S);
    prob  = probSNFE(input)
    V     = solveSNFE(prob,[0.2,1.0,5.0,10.0,20.0])
end
