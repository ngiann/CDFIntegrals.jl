# Following equations taken from Frey and Hinton, 
# Variational Learning in Nonlinear Gaussian Belief Networks,1998
#
# https://www.cs.toronto.edu/~hinton/absps/nlgbn.pdf
#

#---------------------------------------#
# Expectation of ∫ N(x|μ, σ) α⋅Φ(x) dx  #
#---------------------------------------#

M(μ, σ) = Φ(μ/sqrt(1+σ^2))

M(α, μ, σ) = α * M(μ, σ)


#---------------------------------------------------------------------#
# Upper bound to variance of ∫ N(x| μ, σ) ⋅ (α*Φ(x) - M(α, μ, σ))² dx #
#---------------------------------------------------------------------#

V(μ, σ) = Φ(μ/sqrt(1+σ^2)) * (1-Φ(μ/sqrt(1+σ^2))) * (σ^2/(σ^2 + π/2))

V(α, μ, σ) = α * α * V(μ, σ)


#------------------------------------------------------------------------------------#
# Upper bound to expectation of square ∫ N(x|μ, σ) [Φ(x)]² dx = V(μ, σ) + [M(μ, σ)]² #
#------------------------------------------------------------------------------------#

Bslow(μ, σ) = V(μ, σ) + M(μ, σ)^2

Bslow(α, μ, σ) = (α*α) * Bslow(μ, σ)

function B(μ, σ) 
    
    Mμσ = M(μ, σ)

    Mμσ * ( (1-Mμσ) * (σ^2/(σ^2 + π/2)) + Mμσ ) 

end

B(α, μ, σ) = (α*α) * B(μ, σ)