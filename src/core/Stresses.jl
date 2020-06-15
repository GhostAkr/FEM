function calculateStresses(deformations::Array, elasticityMatrix::Matrix)
    stresses = []
    for nodeNum in eachindex(deformations)
        nodeStress = elasticityMatrix * deformations[nodeNum]
        push!(stresses, nodeStress)
    end
    return stresses
end  # stressesKeyword

function calculateVonMises(stresses::Array)
    vonMises = []
    for nodeStress in stresses
        nodeVonMises = sqrt(nodeStress[1]^2 - nodeStress[1] * nodeStress[2] + nodeStress[2]^2 + 3 * nodeStress[3]^2)
        push!(vonMises, nodeVonMises)
    end
    return vonMises
end  # calculateVonMises