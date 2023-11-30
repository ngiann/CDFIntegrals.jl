# Following function is implemented after Wikipedia article
# https://en.wikipedia.org/wiki/List_of_integrals_of_Gaussian_functions
# that presents the solution to the definite integral in (-∞,∞)
# ∫Φ(α + bx)² ϕ(x) dx

#
# This imlements expectation ∫Φ(x)² N(x|μ,σ) dx
#
function J(μ, σ)

    Φ(μ/sqrt(1+σ^2)) - 2*owent(μ/sqrt(1+σ^2), 1/sqrt(1+2σ^2))

end