one_dim_integral(f; lb=-300.0, ub=300.0, dx=1e-4) = mapreduce(f,+,lb:dx:ub)*dx

#-----------------------------------------------------------

function test_integral_M(;μ = μ, σ = σ)

    d = Normal(μ, σ)

    one_dim_integral(x -> pdf(d, x) * Φ(x))

end

#-----------------------------------------------------------

function test_integral_V(;μ = μ, σ = σ)

    d = Normal(μ, σ)

    Mμσ = M(μ, σ)

    one_dim_integral(x ->  pdf(d, x) * (Φ(x) - Mμσ)^2)
    
end

#-----------------------------------------------------------

function test_integral_B(;μ = μ, σ = σ)

    d = Normal(μ, σ)

    one_dim_integral( x-> pdf(d, x) * (Φ(x))^2)
    
end