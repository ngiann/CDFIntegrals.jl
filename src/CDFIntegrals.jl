module CDFIntegrals

    using StatsFuns, HCubature, QuadGK

    include("integral_closed_form.jl")

    
    ϕ(x) = StatsFuns.normpdf(x)
    
    Φ(x) = StatsFuns.normcdf(x)
    
    export M, B, MB
    
end
