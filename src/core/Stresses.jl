function calculateStresses(deformations::Array, elasticityMatrix::Matrix, pars::processPars)
    defsize = size(deformations)
    rowscnt = defsize[1]
    colscnt = defsize[2]
    stresses = zeros(rowscnt, colscnt)
    for node_num in 1:rowscnt
        deformations_vector = [deformations[node_num, value] for value in 1:colscnt]
        node_stress = elasticityMatrix * deformations_vector

        for stressidx in 1:colscnt
            stresses[node_num, stressidx] = node_stress[stressidx]
        end        
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