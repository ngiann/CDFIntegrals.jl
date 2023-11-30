using CDFIntegrals: M, B, V, Φ, ϕ # , Bslow, MB, MBapprox

using Test, Random, HCubature, QuadGK, Printf, Distributions, ProgressMeter

include("test_quadrature.jl")





@testset "CDFIntegrals.jl" begin

    TOL = 1e-7

    numsamples = 30

    rg = MersenneTwister(1) # fix random number generator for reproducinility
   

    #######################################################
    # Verify expectation M(μ, σ) by numerical integration #
    #######################################################

    let 

        greatest_discrepancy = 0.0
        
        @showprogress "Testing M(μ, σ)" for _ in 1:numsamples

            μ, σ = randn(rg)*3, rand(rg)*3

            test = test_integral_M(μ = μ, σ = σ)

            exact = M(μ, σ)

            greatest_discrepancy = max(greatest_discrepancy,  abs(test-exact))
           
        end
        
        @printf("greatest discepancy was %f\n", greatest_discrepancy)
        
        @test greatest_discrepancy < TOL 
        
    end



    ####################################################
    # Verify variance V(μ, σ) by numerical integration #
    ####################################################

    let 
        
        greatest_discrepancy = 0.0
        
        @showprogress "Testing V(μ, σ)" for _ in 1:numsamples

            μ, σ = randn(rg)*3, rand(rg)*3
           
            test = test_integral_V(μ = μ, σ = σ)

            exact = V(μ, σ)

            greatest_discrepancy = max(greatest_discrepancy,  abs(test-exact))
            

        end

        @printf("greatest discepancy was %f\n", greatest_discrepancy)
        
        @test greatest_discrepancy < TOL 

    end


    #################################################################
    # Verify expectation of square B(μ, σ) by numerical integration #
    #################################################################

    let 
        
        greatest_discrepancy = 0.0
        
        @showprogress "Testing B(μ, σ)" for _ in 1:numsamples

            μ, σ = randn(rg)*3, rand(rg)*3
           
            test = test_integral_B(μ = μ, σ = σ)

            exact = B(μ, σ)

            greatest_discrepancy = max(greatest_discrepancy,  abs(test-exact))
            

        end

        @printf("greatest discepancy was %f\n", greatest_discrepancy)
        
        @test greatest_discrepancy < TOL 
        
    end


   
    

    # ##############################################################
    # # Verify that MB(μ, σ) returns M(μ, σ) and B(μ, σ) correctly #
    # ##############################################################

    # let 

    #     @printf("Testing MB(μ, σ)\n")
        
    #     μ, σ = randn(rg)*3, rand(rg)*3
        
    #     m, b = MB(μ, σ)

    #     @test (m == M(μ, σ)) && (b == (B(μ, σ)))
        
    # end


    # ######################################################
    # # Report discrepancy in approximation MBapprox(μ, σ) #
    # ######################################################

    # let 
        
    #     maxdiscrepancy = 0.0

    #     @showprogress "Reporting maximum discrepancy for approximation MBapprox(μ, σ)" for _ in 1:numsamples

    #         μ, σ = randn(rg)*3, rand(rg)*3
        
    #         m, b = MB(μ, σ)
            
    #         mapprox, bapprox = MBapprox(μ, σ)

    #         maxdiscrepancy = max(maxdiscrepancy, abs(mapprox - m), abs(b - bapprox))

    #     end

    #     @printf("\t Maximum discepancy is %f\n", maxdiscrepancy)
        
    # end

end

