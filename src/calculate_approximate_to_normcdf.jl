#
# This script is meant to be used only once in order to establish a good
# and fast approximation to StatsFuns.normcdf via a scaled logistic.
# 
# According to this script the scaling coefficient was determined to be
# equal to 1.7009913570613704 and is used to define the function Φapprox.
#

function calculate_approximate_to_normcdf()

    x = -7:0.1:7

    y = normcdf.(x)

    s(x) = 1/(1+exp(-x))

    s(x, p) = s(x*p)

    objective(p) = sum(abs2.(y - s.(x, p)))

    result = optimize(objective, -10, 10)

    display(result)

    popt = result.minimizer

    x -> s(x, popt), popt

end