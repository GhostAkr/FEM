module multipleIntegral

@enum limits begin
    upper  # Upper limit
    lower  # Lower limit
end

# Multiple integral calculation with cell method. This should be enough to work with squared finite elements.
# f has to be depending on r and s only. So other functions have to be converted to such format before calling this method.
function cellMethod(f:::Function, rLimits::Dict{limits, Number}, sLimits::Dict{limits, Number}, n_r_Sections::Int, n_s_Sections::Int)
    hr = (rLimits.upper - rLimits.lower) / n_r_Sections
    hs = (sLimits.upper - sLimits.lower) / n_s_Sections
    rPoints = collect(StepRange(rLimits.lower, hr, rLimits.upper))
    sPoints = collect(StepRange(sLimits.lower, sr, sLimits.upper))
    print('rPoints: \n', rPoints, '\n')
    print('sPoints: \n', sPoints, '\n')
end

end  # multipleIntegral