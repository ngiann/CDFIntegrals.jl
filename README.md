# CDFIntegrals

[![Build Status](https://github.com/ngiann/CDFIntegrals.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/ngiann/CDFIntegrals.jl/actions/workflows/CI.yml?query=branch%3Amain)

Following equations taken from Frey and Hinton, 
Variational Learning in Nonlinear Gaussian Belief Networks,1998
https://www.cs.toronto.edu/~hinton/absps/nlgbn.pdf

### Expectation
Function `M(α, μ, σ)` returns the exact expectation  of ∫ N(x|μ, σ) α⋅Φ(x) dx where Φ(x) is the standard cumulative normal.

### Variance
Function `V(α, μ, σ)` returns an upper bound to ∫ N(x| μ, σ) ⋅ (α*Φ(x) - M(α, μ, σ))² dx

### Expectation of square
Function `B(α, μ, σ)` returns an upper bound to ∫ N(x|μ, σ) [α⋅Φ(x)]² dx ≤ V(μ, σ) + [M(μ, σ)]² 

