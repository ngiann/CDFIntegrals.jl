using CDFIntegrals
using Test, Random, HCubature, QuadGK, Printf



Bslow(α, μ, σ) = Φ(μ/sqrt(1+σ^2)) * (1-Φ(μ/sqrt(1+σ^2))) * (σ^2/(σ^2 + π/2)) + M(α, μ, σ)^2 


function quadrature_M_1(;μ = μ, σ = σ, α = 1)

    d = Normal(μ, σ)

    integrand(x) = pdf(d, x)*(α * Φ(x))

    hquadrature(integrand, -12, 12)[1]

end


function quadrature_M_2(; μ = μ, σ = σ, α = 1)

    d = Normal(μ, σ)

    integrand(x) = pdf(d, x)*(α * Φ(x))

    quadk(integrand, -12, 12)[1]

end





@testset "CDFIntegrals.jl" begin

    TOL = 1e-3

    numsamples = 10_000

    rg = MersenneTwister(1) # fix random number generator for reproducinility
   

    ##################################################
    # Verify B(μ, σ) with alternative implementation #
    ##################################################

    let 
        
        FLAG = true

        for _ in 1:numsamples
            
            α, μ, σ = rand()*3, randn(rg)*3, rand(rg)*3
            
            FLAG *= (abs(Bslow(α, μ, σ) - B(α, μ, σ)) < TOL)

        end
        
        @test FLAG

    end


    ###########################################
    # Verify M(μ, σ) by numerical integration #
    ###########################################

    let 
        @printf("Testing expectation of Φ(x)\n")

        test1 = test_quadrature_M_1(μ = μ, σ = σ, α = α)
        test2 = test_quadrature_M_2(μ = μ, σ = σ, α = α)
 
        M(α, μ, σ)
           
        if abs(test1-test2)>TOL || abs(test2-test3)>TOL || abs(test1-test3) > TOL
            @printf("test1 = %f, test2 = %f, test3 = %f\n", test1, test2, test3)
            error("Tolerance exceeded\n")    
        end
    end

end

