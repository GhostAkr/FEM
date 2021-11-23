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
    von_mises = []
    stresses_size = size(stresses)
    stresses_rows = stresses_size[1]
    stresses_cols = stresses_size[2]

    for stress_node_idx in 1:stresses_rows

        node_von_mises = 0
        if stresses_cols == 3  # 2D case
            sxx = stresses[stress_node_idx, 1]
            syy = stresses[stress_node_idx, 2]
            sxy = stresses[stress_node_idx, 3]

            node_von_mises = sqrt(sxx^2 - sxx * syy + syy^2 + 3 * sxy^2)
        elseif stresses_cols == 6  # 3D case
            sxx = stresses[stress_node_idx, 1]
            syy = stresses[stress_node_idx, 2]
            szz = stresses[stress_node_idx, 3]
            sxy = stresses[stress_node_idx, 4]
            syz = stresses[stress_node_idx, 5]
            szx = stresses[stress_node_idx, 6]

            node_von_mises = sqrt(0.5 * ((sxx - syy)^2 + (syy - szz)^2 + (szz - sxx)^2 + 
                6 * (syz^2 + szx^2 + sxy^2)))
        else
            @error("Error while calculating von Mises stresses in calculateVonMises()")
            return nothing
        end

        push!(von_mises, node_von_mises)
    end
    return von_mises
end  # calculateVonMises