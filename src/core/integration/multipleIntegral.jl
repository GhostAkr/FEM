module multipleIntegral

include("multipleIntegralMethods.jl")

function integrateMatrix(matrix::Array{Number, 2}, method::integrateMethods, rLim::limits, sLim::limits)
    resultMatrix = Array{typeof(matrix[1, 1]), 2}(undef, size(matrix, 1), size(matrix, 2))
    if method == cell begin
        M = 3200
        N = 3200
        for i in size(matrix, 1)
            for j in size(matrix, 2)
                resultMatrix[i, j] = cellMethod(matrix[i, j], rLim, sLim, M, N)
            end
        end
    end
        
end

end  # multipleIntegral