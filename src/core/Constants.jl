# TODO: Make it for different problems (e.g. plane stress etc.)
"""
    elasticityMatrix(youngMod, poissRatio)

Calculates plane strain elasticity matrix with given Young's modulus and Poisson's ratio.

# Arguments
- `youngMod::Float64`: Young's modulus;
- `poissRatio::Float64`: Poisson's ratio.
"""
function elasticityMatrix(youngMod, poissRatio, type::Int)
    nu = poissRatio
    E = youngMod
    if type == 1
        resMatr = [1    nu / (1 - nu)   0
                    nu / (1 - nu)   1   0
                    0   0   (1 - 2 * nu) / (2 * (1 - nu))]
        resMatr *= (E * (1 - nu) / ((1 + nu) * (1 - 2 * nu)))
    elseif type == 2
        resMatr = [1 nu 0
                    nu 1 0
                    0 0 (1 - nu) / 2]
        resMatr *= (E / (1 - nu^2))
    else
        println("Given model type in elasticityMatrix is not supported")
        return nothing
    end
    return resMatr
end  # elasticityMatrix