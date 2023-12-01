# Following equations taken from Frey and Hinton, 
# Variational Learning in Nonlinear Gaussian Belief Networks,1998
#
# https://www.cs.toronto.edu/~hinton/absps/nlgbn.pdf
#

#-----------------------------------[1]-#
# Expectation of ∫ N(x|μ, σ) Φ(x) dx    #
#---------------------------------------#

M(μ, σ) = Φ(μ/sqrt(1+σ^2))


#-----------------------------------[2]-#
# Expectation of ∫ N(x|μ, σ) Φ(x*s) dx  #
#---------------------------------------#

M(μ, σ, s) = M(μ*s, σ*s)


#---------------------------------------------[3]-#
# Expectation of ∫ N(x|μ, σ) (Φ(x) - M(μ,σ))² dx  #
#-------------------------------------------------#

V(μ, σ) = J(μ, σ) - M(μ,σ)^2


#-----------------------------------[4]-#
# Expectation of ∫ N(x|μ, σ) [Φ(x)]² dx #
#---------------------------------------#

B(μ, σ) = J(μ, σ)


#-------------------------------------[5]-#
# Expectation of ∫ N(x|μ, σ) [Φ(x*s)]² dx #
#-----------------------------------------#

B(μ, σ, s) = B(μ*s, σ*s)


#-----------------------------------[6]-#
# Expectation of ∫ N(x|μ, σ) M(x, s) dx #
#                                       #
# where M(x,s) = Φ(x/sqrt(1+s^2))       #
#---------------------------------------#

function MM(μ, σ, s)
   
    scale = sqrt(1+s^2)

    M(μ, σ,  1/scale)

end


#--------------------------------------[7]-#
# Expectation of                           # 
# ∫ N(x|μ, σ) [M(x, s)]² dx =              #
# ∫ N(x|μ, σ) [Φ(x/sqrt(1+s^2))]² dx       #
#------------------------------------------#

function MM²(μ, σ, s)
   
    scale = sqrt(1+s^2)

    B(μ, σ, 1/scale)

end



#-------------------------------------[8]-#
# Expectation of ∫ N(x|μ, Σ) M(vᵀx, s) dx #
#                                         #
# where M(x,s) = Φ(x/sqrt(1+s^2))         #
#-----------------------------------------#

function MM(μ, Σ, s, v)
   
    scale = sqrt(1+s^2)

    M(dot(v, μ), sqrt(dot(v, Σ*v)), 1/scale) # consult Barber, Eq. (28.5.8)

end



#----------------------------------------------------------------#
# Upper bound to variance of ∫ N(x| μ, σ) ⋅ (Φ(x) - M(μ, σ))² dx #
#----------------------------------------------------------------#

function Vup(μ, σ) 
    
    Mμσ = M(μ, σ)

    Mμσ * (1-Mμσ) * (σ^2/(σ^2 + π/2)) # note: part of the expression in Bup(μ, σ)

end


#-------------------------------------------------------------#
# Upper bound to expectation of square ∫ N(x|μ, σ) [Φ(x)]² dx #
#-------------------------------------------------------------#

function Bup(μ, σ) 
    
    Mμσ = M(μ, σ)

    Mμσ * ( (1-Mμσ) * (σ^2/(σ^2 + π/2)) + Mμσ ) 

end



# #----------------------------------------------------------------------#
# # Convenience function that returns both M(μ, σ) and B(μ, σ) in one go #
# #----------------------------------------------------------------------#

# function MB(μ, σ) 
    
#     σ² = σ^2

#     Mμσ = Φ(μ/sqrt(1+σ²))

#     Bμσ = Mμσ * ( (1-Mμσ) * (σ²/(σ² + π/2)) + Mμσ ) 

#     Mμσ, Bμσ

# end


# #----------------------------------------------------------------------------#
# # Approximate, faster version of Φ and MB                                    #
# # Constant 1.7009913570613704 obtained from calculate_approximate_to_normcdf #
# #----------------------------------------------------------------------------#

# # sigm(x) = 1.0 / (1.0 + exp(-x))

# # Φapprox(x) = sigm(x * 1.7009913570613704) 

# Φapprox(x) = 1.0 / (1.0 + exp(-x*1.7009913570613704)) 

# function MBapprox(μ, σ) 
    
#     σ² = σ^2

#     Mμσ = Φapprox(μ/sqrt(1+σ²))

#     Bμσ = Mμσ * ( (1-Mμσ) * (σ²/(σ² + π/2)) + Mμσ ) 

#     Mμσ, Bμσ

# end