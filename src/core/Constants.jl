# TODO: Make it for different problems (e.g. plane stress etc.)
"""
    elasticityMatrix(youngMod, poissRatio)

Calculates plane strain elasticity matrix with given Young's modulus and Poisson's ratio.

# Arguments
- `youngMod::Float64`: Young's modulus;
- `poissRatio::Float64`: Poisson's ratio.
"""
function elasticityMatrix(youngMod, poissRatio)
    nu = poissRatio
    E = youngMod
    resMatr = [1    nu / (1 - nu)   0
                nu / (1 - nu)   1   0
                0   0   (1 - 2 * nu) / (2 * (1 - nu))]
    resMatr *= (E * (1 - nu) / ((1 + nu) * (1 - 2 * nu)))
end  # elasticityMatrix