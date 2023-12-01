module CDFIntegrals

    using StatsFuns, HCubature, QuadGK, IrrationalConstants, SpecialFunctions, LinearAlgebra

    include("integral_closed_form.jl")

    include("owent.jl") # code written by Andy Gough, uploaded in julia forum

    include("EPhiSquared.jl")
    
    ϕ(x) = StatsFuns.normpdf(x)
    
    Φ(x) = StatsFuns.normcdf(x)
    
    export Φ

    export M, B, V

    export Bup
    
end
