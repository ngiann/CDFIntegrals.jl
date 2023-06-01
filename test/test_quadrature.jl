one_dim_integral(f; lb=-100.0, ub=100.0, dx=1e-4) = mapreduce(f,+,lb:dx:ub)*dx

#-----------------------------------------------------------

function test_integral_M(;μ = μ, σ = σ, α = α)

    d = Normal(μ, σ)

    one_dim_integral(x -> pdf(d, x) * α * Φ(x))

end

#-----------------------------------------------------------

function test_integral_V(;μ = μ, σ = σ, α = α)

    d = Normal(μ, σ)

    one_dim_integral(x ->  pdf(d, x) * (α * Φ(x) - M(α, μ, σ))^2)
    
end

#-----------------------------------------------------------

function test_integral_B(;μ = μ, σ = σ, α = α)

    d = Normal(μ, σ)

    one_dim_integral( x-> pdf(d, x) * (α * Φ(x))^2)
    
end