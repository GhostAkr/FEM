function calculateStresses(deformations::Array, elasticityMatrix::Matrix, pars::processPars)
    stresses = []
    for nodeNum in 1:size(deformations)[1]
        deformationsVector = [deformations[nodeNum, value] for value in 1:size(deformations)[2]]
        nodeStress = elasticityMatrix * deformationsVector
        # append!(nodeStress, pars.materialProperties[poisC] * (nodeStress[1] + nodeStress[2]))
        push!(stresses, nodeStress)
    end
    return stresses
end  # stressesKeyword

function calculateVonMises(stresses::Array)
    vonMises = []
    for nodeStress in stresses
        nodeVonMises = sqrt(nodeStress[1]^2 - nodeStress[1] * nodeStress[2] + nodeStress[2]^2 + 3 * nodeStress[3]^2)
        sxx = nodeStress[1]  # Stress by XX
        syy = nodeStress[2]  # Stress by YY
        sxy = nodeStress[3]  # Stress by XY
        # szz = nodeStress[4]  # Stress by ZZ
        # nodeVonMises = sqrt(0.5 * ((sxx - syy)^2 + (syy - szz)^2 + (szz - sxx)^2 + 6 * sxy^2))
        push!(vonMises, nodeVonMises)
    end
    return vonMises
end  # calculateVonMises