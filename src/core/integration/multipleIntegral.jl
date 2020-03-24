module multipleIntegral

export upper, lower

@enum limits begin
    upper  # Upper limit
    lower  # Lower limit
end

# Multiple integral calculation with cell method. This should be enough to work with squared finite elements.
# f has to be depending on r and s only. So other functions have to be converted to such format before calling this method.
# This method is very simple but not effective. So need to use another one to improve speed.
function cellMethod(f, rLimits::Dict{limits, Real}, sLimits::Dict{limits, Real}, n_r_Sections::Int, n_s_Sections::Int)
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

end  # multipleIntegral