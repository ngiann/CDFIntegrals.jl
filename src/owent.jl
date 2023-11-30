# Owen's T Function
#    Written by Andy Gough; August 2021 
#    Rev 1.09
#    MIT License 
#
# dependencies
#   StatsFuns
#   IrrationalConstants
#   SpecialFunctions
#   LinearAlgebra
#
# HISTORY 
# In the past 20 or so years, most implementations of Owen's T function have followed the algorithms given in "Fast and accurate Calculation of Owen's 
# T-Function", by M. Patefield and D. Tandy, Journal of Statistical Software, 5 (5), 1 - 25 (2000)
#
# Six algorithms were given, and which is was used depends on the values of (h,a) 
#
# T1: first m terms of series expansion of Owen (1956)
# T2: approximates 1/(1+x^2) by power series expansion up to order 2m
# T3: approximates 1/(1+x^2) by chebyshev polynomials of degree 2m in x
# T4: new expression for zi from T2
# T5: Gauss 2m-point quadrature; 30 figures accuracy with m=48 (p. 18)
# T6: For when a is very close to 1, use formula derived from T(h,1) = 1/2 Φ(h)[1-Φ(h)]
#
# They developed code for these algorithms on a DEC VAX 750. The VAX 750 came out in 1980 and had a processor clock speed of 3.125 MHz.
#
# The reason for 6 algorithms was to speed up the function when possible, with T1 being faster than T2, T2 faster than T3, etc.
#
# THIS FUNCTION 
#   A native Julia implementation, based on the equations in the paper. The FORTRAN source code was not anlayzed, translated, or used. This is a new
# implementation that takes advantages of Julia's unique capabilities (and those of modern computers). 
#
# T1 through T4 are not implemented. Instead, if a < 0.999999, T5 is used to calculate Owen's T (using 48 point Gauss-Legendre quadrature)
# For 0.999999 < a < 1.0, T6 is implemented. 
#
# REFERENCES
# [1] "Fast and accurate Calculation of Owen's T-Function", by M. Patefield and D. Tandy, Journal of Statistical Software, 5 (5), 1 - 25 (2000)
# [2] "Tables for Computing Bivariate Normal Probabilities", by Donald P. Owen, The Annals of Mathematical Statistics, Vol. 27, No. 4 (Dec 1956), pp. 1075-1090
# 
# Partial Derivatives (FYI)
# D[owent[x,a],x] = -exp(-0.5*x^2)*erf(a*x/sqrt2)/(2*sqrt2π)
# D[owent[x,a],a] = exp(-0.5*(1+a^2)*(x^2))/((1+a^2)*2π)
#
@doc raw"""
owent(h,a) : Returns the value of Owen's T function for (h,a)

Owen's T function:
```math
T(h,a) = \frac{1}{2\pi } \int_{0}^{a} \frac{e^{-\frac{1}{2}h^2(1+x^2)}}{1+x^2}dx\quad(-\infty < h,a < +\infty)
```

For *h* and *a* > 0, *T(h,a)* gives the volume of the uncorrelated bivariate normal distribution with zero mean and unit variance 
over the area from *y = ax* and *y = 0* and to the right of *x = h*. 

EXAMPLE:
```
julia> owent(0.0625, 0.025)
0.003970281304296922
```

Worst case accuracy is about 2E-16. 
"""
function owent(h::T1,a::T1) where {T1 <: Real}

    #*********************
    # shortcut evaluations 
    #*********************

    if h < zero(h) 
        return owent(abs(h),a)
    end 

    if h==zero(h) 
        return T1(atan(a)*inv2π)
    end

    if a < zero(a) 
        return -owent(h,abs(a))
    end 

    if a==zero(a)
        return zero(a)
    end 

    if a==one(a) 
        return T1(0.125*erfc(-h*invsqrt2)*erfc(h*invsqrt2))
    end 

    if a==Inf 
        return T1(0.25*erfc(sqrt(h^2)*invsqrt2))
    end 

    # below reduces the range from -inf < h,a < +inf to h ≥ 0, 0 ≤ a ≤ 1
    if a > one(a) 
        return 0.25*(erfc(-h*invsqrt2) + erfc(-a*h*invsqrt2)) - 0.25*erfc(-h*invsqrt2)*erfc(-a*h*invsqrt2) - owent(a*h,one(a)/a)
    end

    # calculate Owen's T 

    if a ≤ 0.999999 

        t2(h,a,x) = inv4π*a*exp(-0.5*(h^2)*(one(h)+(a*x)^2))/(one(h)+(a*x)^2)

        x = [-0.9987710072524261, -0.9935301722663508, -0.9841245837228269, -0.9705915925462473, -0.9529877031604309, -0.9313866907065543, -0.9058791367155696
        , -0.8765720202742479, -0.8435882616243935, -0.8070662040294426, -0.7671590325157404, -0.7240341309238146, -0.6778723796326639, -0.6288673967765136
        , -0.5772247260839727, -0.5231609747222331, -0.4669029047509584, -0.4086864819907167, -0.34875588629216075, -0.28736248735545555, -0.22476379039468905
        , -0.1612223560688917, -0.0970046992094627, -0.03238017096286937, 0.03238017096286937, 0.0970046992094627, 0.1612223560688917, 0.22476379039468905
        , 0.28736248735545555, 0.34875588629216075, 0.4086864819907167, 0.4669029047509584, 0.5231609747222331, 0.5772247260839727, 0.6288673967765136
        , 0.6778723796326639, 0.7240341309238146, 0.7671590325157404, 0.8070662040294426, 0.8435882616243935, 0.8765720202742479, 0.9058791367155696
        , 0.9313866907065543, 0.9529877031604309, 0.9705915925462473, 0.9841245837228269, 0.9935301722663508, 0.9987710072524261]

        w = [0.0031533460523059122, 0.0073275539012762885, 0.011477234579234613, 0.015579315722943824, 0.01961616045735561, 0.023570760839324363
        , 0.027426509708356944, 0.031167227832798003, 0.03477722256477052, 0.038241351065830737, 0.04154508294346467, 0.0446745608566943, 0.04761665849249045
        , 0.05035903555385445, 0.05289018948519363, 0.055199503699984116, 0.05727729210040322, 0.05911483969839564, 0.06070443916589387, 0.06203942315989268
        , 0.06311419228625402, 0.06392423858464813, 0.06446616443594998, 0.06473769681268386, 0.06473769681268386, 0.06446616443594998, 0.06392423858464813
        , 0.06311419228625402, 0.06203942315989268, 0.06070443916589387, 0.05911483969839564, 0.05727729210040322, 0.055199503699984116
        , 0.05289018948519363, 0.05035903555385445, 0.04761665849249045, 0.0446745608566943, 0.04154508294346467, 0.038241351065830737, 0.03477722256477052
        , 0.031167227832798003, 0.027426509708356944, 0.023570760839324363, 0.01961616045735561, 0.015579315722943824, 0.011477234579234613, 0.0073275539012762885
        , 0.0031533460523059122]
        
        towen = dot(w, t2.(h,a,x))

        return T1(towen)
    else
        # a > 0.999999, T6 from paper (quadrature using QuadGK would also work, but be slower)

        j = 0.5*erfc(-h*invsqrt2)
        k = atan((one(a)-a)/(one(a)+a))
        towen = 0.5*j*(one(h)-j)-inv2π*k*exp((-0.5*(one(a)-a)*h^2)/k)

        return T1(towen)
    end 
end 

owent(h::Real, a::Real) = owent(promote(h,a)...)