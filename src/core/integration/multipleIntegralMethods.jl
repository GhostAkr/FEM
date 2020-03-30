export upper, lower, cell, limits

@enum limits begin
    upper  # Upper limit
    lower  # Lower limit
end

@enum integrateMethods begin
    cell  # Cell integration method
end

# Multiple integral calculation with cell method. This should be enough to work with squared finite elements.
# f has to be depending on r and s only. So other functions have to be converted to such format before calling this method.
# This method is very simple but not effective. So need to use another one to improve speed.
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

#Integrating matrix of functions
# F should return matrix depending on two variables, See StiffnessMatrix.jl for details
function cellMethodMatrix(F::Function, rLimits::Dict{limits, Real}, sLimits::Dict{limits, Real}, n_r_Sections::Int, n_s_Sections::Int)
    hr = (rLimits[upper] - rLimits[lower]) / n_r_Sections
    hs = (sLimits[upper] - sLimits[lower]) / n_s_Sections
    rPoints = collect(Float64, rLimits[lower]:hr:rLimits[upper])
    sPoints = collect(Float64, sLimits[lower]:hs:sLimits[upper])
    sum = 0
    nOfRows = size(F(1, 1))[1]
    nOfCols = size(F(1, 1))[2]
    resultMatrix = Matrix{Float64}(undef, nOfRows, nOfCols)
    # Going through all matrix elements. Integrating each one.
    for k in 1:nOfRows
        for l in 1:nOfCols
            # Integrating [k, l] matrix element
            for i in rPoints
                #println("i = ", i)
                buffSum = 0
                for j in sPoints
                    buffSum += (F(i, j)[k, l] * hr * hs)
                    #print("Element #", i * nOfCols + j)
                end
                sum += buffSum
            end
            resultMatrix[k, l] = sum
            #println("Matrix element ", k, ", ", l)
        end
    end
    return resultMatrix
end