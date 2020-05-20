export upper, lower, cell, limits

"""
    limits

Possible intergral limits:
1. `upper`: upper limit;
2. `lower`: lower limit.
"""
@enum limits begin
    upper
    lower
end

"""
    cellMethod(f::Function, rLimits::Dict{limits, Real}, sLimits::Dict{limits, Real}, n_r_Sections::Int, n_s_Sections::Int)

Multiple integral calculation with cell method. f has to be depended on r and s only. 
So other functions have to be converted to such format before calling this method.
This method is very simple but not effective. This method was implemented for tests, so it's recommended to use Gauss integration instead.

# Arguments
- `f`: function that needs to be integrated, need to be depended on 2 variables: ``f = f(r, s)``;
- `rLimits::Dict{limits, Real}`: limits for the first variable;
- `sLimits::Dict{limits, Real}`: limits for the second variable;
- `n_r_Sections::Int`: number of sections on ``r`` axis;
- `n_s_Sections::Int`: number of sections on ``s`` axis.
"""
function cellMethod(f::Function, rLimits::Dict{limits, Real}, sLimits::Dict{limits, Real}, n_r_Sections::Int, n_s_Sections::Int)
    hr = (rLimits[upper] - rLimits[lower]) / n_r_Sections
    hs = (sLimits[upper] - sLimits[lower]) / n_s_Sections
    rPoints = collect(Float64, rLimits[lower]:hr:rLimits[upper])
    sPoints = collect(Float64, sLimits[lower]:hs:sLimits[upper])
    sum = 0
    for i in rPoints
        buffSum = 0
        for j in sPoints
            buffSum += (f(i, j) * hr * hs)
        end
        sum += buffSum
    end
    return sum
end

# Integrating matrix of functions
# F should return matrix depending on two variables, See StiffnessMatrix.jl for details
# !IMPORTANT! Probably incorrect. For now gauss method can be used (and it's definitely more efficient).
# function cellMethodMatrix(F::Function, rLimits::Dict{limits, Real}, sLimits::Dict{limits, Real}, n_r_Sections::Int, n_s_Sections::Int)
#     hr = (rLimits[upper] - rLimits[lower]) / n_r_Sections
#     hs = (sLimits[upper] - sLimits[lower]) / n_s_Sections
#     rPoints = collect(Float64, rLimits[lower]:hr:rLimits[upper])
#     sPoints = collect(Float64, sLimits[lower]:hs:sLimits[upper])
#     sum = 0
#     nOfRows = size(F(1, 1))[1]
#     nOfCols = size(F(1, 1))[2]
#     resultMatrix = Matrix{Float64}(undef, nOfRows, nOfCols)
#     # Going through all matrix elements. Integrating each one.
#     for k in 1:nOfRows
#         for l in 1:nOfCols
#             # Integrating [k, l] matrix element
#             for i in rPoints
#                 #println("i = ", i)
#                 buffSum = 0
#                 for j in sPoints
#                     buffSum += (F(i, j)[k, l] * hr * hs)
#                     #print("Element #", i * nOfCols + j)
#                 end
#                 sum += buffSum
#             end
#             resultMatrix[k, l] = sum
#             #println("Matrix element ", k, ", ", l)
#         end
#     end
#     return resultMatrix
# end

"""
    gaussMethod(F::Function, intOrder::Int)

Multiple integral calculation with Gauss method. F has to be depended on r and s only. 
So other functions have to be converted to such format before calling this method.
Interval of integration should be equal to [-1; 1].
Supported integration orders: 2.

# Arguments
- `F::Function`: function that needs to be integrated, need to be depended on 2 variables: ``F = F(r, s)``;
- `intOrder::Int`: integration order.
"""
function gaussMethod(F::Function, intOrder::Int)
    rArray = Array  # Array of integration points by r coordinate
    sArray = Array  # Array of integration points by s coordinate
    weights = Array  # Array of integration weights
    if intOrder == 2
        r = [-1 / sqrt(3), 1 / sqrt(3)]
        s = [-1 / sqrt(3), 1 / sqrt(3)]
        weights = [1, 1]
    else
        println("That integration order is not supported")
    end
    resultSum = 0
    for i in 1:intOrder
        for j in 1:intOrder
            resultSum += (weights[i] * weights[j] * F(r[i], s[j]))
        end
    end
    return resultSum
end  # gaussMethod

"""
    gaussMethodMatrix(F::Function, intOrder::Int)

Integrate matrix of functions depending on 2 variables with Gauss method.

# Arguments
- `F::Function`: functions returning matrix of functions depending on 2 variables;
- `intOrder::Int`: integration order.
"""
function gaussMethodMatrix(F::Function, intOrder::Int)
    rArray = Array  # Array of integration points by r coordinate
    sArray = Array  # Array of integration points by s coordinate
    weights = Array  # Array of integration weights
    if intOrder == 2
        r = [-1 / sqrt(3), 1 / sqrt(3)]
        s = [-1 / sqrt(3), 1 / sqrt(3)]
        weights = [1, 1]
    else
        println("That integration order is not supported")
    end
    nOfRows = size(F(1, 1))[1]
    if length(size(F(1, 1))) == 1
        nOfCols = 1
    else
        nOfCols = size(F(1, 1))[2]
    end
    resultMatrix = Matrix{Float64}(undef, nOfRows, nOfCols)
    for k in 1:nOfRows
        for l in 1:nOfCols
            # Integrating [k, l] element of source matrix
            resultSum = 0
            for i in 1:intOrder
                for j in 1:intOrder
                    resultSum += (weights[i] * weights[j] * F(r[i], s[j])[k, l])
                end
            end
            resultMatrix[k, l] = resultSum
        end
    end
    return resultMatrix
end  # gaussMethodMatrix

"""
    gauss1DMethodMatrix(F::Function, intOrder::Int)

Integrate matrix of functions depending on 1 variable with Gauss method.

# Arguments
- `F::Function`: functions returning matrix of functions depending on 1 variable;
- `intOrder::Int`: integration order.
"""
function gauss1DMethodMatrix(F::Function, intOrder::Int)
    if intOrder == 2
        x = [-1 / sqrt(3), 1 / sqrt(3)]
        weights = [1, 1]
    else
        println("That integration order is not supported")
    end
    nOfRows = size(F(1))[1]
    if length(size(F(1))) == 1
        nOfCols = 1
    else
        nOfCols = size(F(1))[2]
    end
    resultMatrix = Matrix{Float64}(undef, nOfRows, nOfCols)
    for k in 1:nOfRows
        for l in 1:nOfCols
            # Integrating [k, l] element of source matrix
            resultSum = 0
            for i in 1:intOrder
                resultSum += (weights[i] * F(x[i])[k, l])
            end
            resultMatrix[k, l] = resultSum
        end
    end
    return resultMatrix
end  # gauss1DMethodMatrix