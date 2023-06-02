# CDFIntegrals

[![Build Status](https://github.com/ngiann/CDFIntegrals.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/ngiann/CDFIntegrals.jl/actions/workflows/CI.yml?query=branch%3Amain)

Following equations taken from Frey and Hinton, 
Variational Learning in Nonlinear Gaussian Belief Networks,1998
https://www.cs.toronto.edu/~hinton/absps/nlgbn.pdf

### Expectation
Exported function `M(α, μ, σ)` returns the exact expectation  of ∫ N(x|μ, σ) α⋅Φ(x) dx where Φ(x) is the standard cumulative normal.

### Variance
Non-exported function `V(α, μ, σ)` returns an upper bound to ∫ N(x| μ, σ) ⋅ (α*Φ(x) - M(α, μ, σ))² dx.

### Expectation of square
Exported function `B(α, μ, σ)` returns an upper bound to ∫ N(x|μ, σ) [α⋅Φ(x)]² dx ≤ V(μ, σ) + [M(μ, σ)]² = B(α, μ, σ).

### Example

```
using CDFIntegrals
using PyPlot # must be indepedently installed

plot(x, M.(x,σ),"blue",label="mean of noisy cdfs")
plot(x, M.(x,σ) .+ CDFIntegrals.V.(x,σ),"--r",label="variance")
plot(x, M.(x,σ) .- CDFIntegrals.V.(x,σ),"--r")
plot(x, CDFIntegrals.Φ.(x),"g",label="cdf")

legend()
```
![example](cdf_with_interval.png)