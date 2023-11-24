using CDFIntegrals: M, B, V, Φ, ϕ, Bslow, MB, MBapprox

using Test, Random, HCubature, QuadGK, Printf, Distributions, ProgressMeter

include("test_quadrature.jl")





@testset "CDFIntegrals.jl" begin

    TOL = 1e-8

    numsamples = 100

    rg = MersenneTwister(1) # fix random number generator for reproducinility
   

    ##########################################################
    # Verify expectation M(α, μ, σ) by numerical integration #
    ##########################################################

    let 

        greatest_discrepancy = 0.0
        
        @showprogress "Testing M(α, μ, σ)" for _ in 1:numsamples

            α, μ, σ = rand(rg)*3, randn(rg)*3, rand(rg)*3

            test = test_integral_M(μ = μ, σ = σ, α = α)

            exact = M(α, μ, σ)

            greatest_discrepancy = max(greatest_discrepancy,  abs(test-exact))
           
        end
        
        @printf("greatest discepancy was %f\n", greatest_discrepancy)
        
        @test greatest_discrepancy < TOL 
        
    end



    ##########################################################
    # Verify upper bound V(α, μ, σ) by numerical integration #
    ##########################################################

    let 
        
        is_upper_bound = true
        
        @showprogress "Checking that V(α, μ, σ) is upper bound" for _ in 1:numsamples

            α, μ, σ = rand(rg)*3, randn(rg)*3, rand(rg)*3
           
            test = test_integral_V(μ = μ, σ = σ, α = α)

            exact = V(α, μ, σ)

            is_upper_bound *= (exact >= test)
            

        end

        @test is_upper_bound
    end



    #####################################################
    # Verify B(α, μ, σ) with alternative implementation #
    #####################################################

    let 

        @printf("Testing B(α, μ, σ) with alternative implementation... ")
        
        greatest_discrepancy = 0.0

        for _ in 1:numsamples
            
            α, μ, σ = rand(rg)*3, randn(rg)*3, rand(rg)*3
            
            greatest_discrepancy = max(greatest_discrepancy, abs(Bslow(α, μ, σ) - B(α, μ, σ)))
            greatest_discrepancy = max(greatest_discrepancy, abs(α * Bslow(μ, σ) - α * B(μ, σ)))

        end
        
        @printf("greatest discrepancy was %f\n", greatest_discrepancy)

        @test greatest_discrepancy < TOL

    end


    ###############################################################
    # Verify B(α, μ, σ) is an upper bound to expectation of square
    ###############################################################

    let 
        
        is_upper_bound = true
        
        @showprogress "Checking that B(α, μ, σ) is upper bound" for _ in 1:numsamples

            α, μ, σ = rand(rg)*3, randn(rg)*3, rand(rg)*3
           
            test = test_integral_B(μ = μ, σ = σ, α = α)

            exact = B(α, μ, σ)

            is_upper_bound *= (exact + TOL >= test)
            

        end

        @test is_upper_bound == true
        
    end


    ###############################################################
    # Verify that MB(μ, σ) returns M(μ, σ) and B(μ, σ) correctly
    ###############################################################

    let 

        @printf("Testing MB(μ, σ)\n")
        
        μ, σ = randn(rg)*3, rand(rg)*3
        
        m, b = MB(μ, σ)

        @test (m == M(μ, σ)) && (b == (B(μ, σ)))
        
    end


    ###############################################################
    # Report discrepancy in approximation MBapprox(μ, σ)
    ###############################################################

    let 
        
        maxdiscrepancy = 0.0

        @showprogress "Reporting maximum discrepancy for approximation MBapprox(μ, σ)" for _ in 1:numsamples

            μ, σ = randn(rg)*3, rand(rg)*3
        
            m, b = MB(μ, σ)
            
            mapprox, bapprox = MBapprox(μ, σ)

            maxdiscrepancy = max(maxdiscrepancy, abs(mapprox - m), abs(b - bapprox))

        end

        @printf("\t Maximum discepancy is %f\n", maxdiscrepancy)
        
    end

end

