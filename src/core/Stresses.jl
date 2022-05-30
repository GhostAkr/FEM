using Quad4Pts
using Quad8Pts
using Iso8Pts3D
using MultipleIntegral

"""
    calculateStresses(deformations::Array, elasticityMatrix::Matrix)

Calculate stresses for usual elastic mechanical problem.

Param `deformations::Array`:
    Elements of this array are arrays of deformations in different directions for given 
    node. Length of each such array depends on model. For example in 2D model we have only 3
    deformation directions and in 3D model we have 6 of these.

# Arguments
- `deformations::Array`: array of deformations;
- `elasticityMatrix::Matrix`: elasticity matrix.
"""
function calculateStresses(deformations::Array, elasticityMatrix::Matrix)
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
end  # calculateStresses

"""
    nonloc_stresses_integrand(r_source, s_source, t_source, r_impact, s_impact, 
        t_impact, x_source::Array{Float64}, y_source::Array{Float64}, 
        z_source::Array{Float64}, x_impact::Array{Float64}, y_impact::Array{Float64}, 
        z_impact::Array{Float64}, displacements::Array, impactdist::Number, 
        elasticitymatrix::AbstractArray, elemtype::FiniteElement)

Computes integrand to calculate non-local part of stresses. `displacements` here represents
displacements of nodes in impact element only.

# Arguments
- `r_source`: r-coordinate of node in source element;
- `s_source`: s-coordinate of node in source element;
- `t_source`: t-coordinate of node in source element;
- `r_impact`: r-coordinate of node in impact element;
- `s_impact`: s-coordinate of node in impact element;
- `t_impact`: t-coordinate of node in impact element;
- `x_source::Array{Float64}`: x-coordinates of nodes in source element;
- `y_source::Array{Float64}`: y-coordinates of nodes in source element;
- `z_source::Array{Float64}`: z-coordinates of nodes in source element;
- `x_impact::Array{Float64}`: x-coordinates of nodes in impact element;
- `y_impact::Array{Float64}`: y-coordinates of nodes in impact element;
- `z_impact::Array{Float64}`: z-coordinates of nodes in impact element;
- `displacements::Array`: displacements vector;
- `impactdist::Number`: impact distance;
- `elasticitymatrix::AbstractArray`: elasticity matrix;
- `elemtype::FiniteElement`: type of finite element.
"""
function nonloc_stresses_integrand(r_source, s_source, t_source, r_impact, s_impact, 
    t_impact, x_source::Array{Float64}, y_source::Array{Float64}, z_source::Array{Float64}, 
    x_impact::Array{Float64}, y_impact::Array{Float64}, z_impact::Array{Float64},
    displacements::Array, impactdist::Number, elasticitymatrix::AbstractArray, 
    elemtype::FiniteElement
)
    gauss_node_deformations = gradMatr(r_impact, s_impact, t_impact, 
        x_impact, y_impact, z_impact, elemtype) * displacements

    integrvec = elasticitymatrix * gauss_node_deformations

    # Calculating impact function
    sourcept_glob = conv_loc_to_glob(r_source, s_source, t_source, x_source, y_source, 
        z_source, elemtype)
    impactpt_glob = conv_loc_to_glob(r_impact, s_impact, t_impact, x_impact, y_impact, 
        z_impact, elemtype)
    impactvec = vecfrompoints_3d(sourcept_glob, impactpt_glob)
    currdist = veclength_3d(impactvec)
    normfact = 1 / (8 * pi * impactdist^3)
    impact = nonloc_gaussimpact(normfact, impactdist, currdist)

    integrvec .*= impact

    return integrvec
end

"""
    calulate_stresses_3d_nonloc_depr(deformations::Array, displacements::Array, 
        elasticitymatrix::Matrix, beta_loc::Number, beta_nonloc::Number, neighbours::Array, 
        impactdist::Number, pars::ProcessPars, intorder::Int, elemtype::FiniteElement)

Calculate stresses in non-local case taking neighbours impact into account.

# Arguments
- `deformations::Array`: deformations in nodes for each element;
- `displacements::Array`: displacements vector;
- `elasticitymatrix::Matrix`: elasticity matrix;
- `beta_loc::Number`: coefficient which defines impact of local part of stiffness matrix;
- `beta_nonloc::Number`: coefficient which defines impact of non-local part of stiffness 
    matrix.
- `neighbours::Array`: array of neighbours for each element;
- `impactdist::Number`: impact distance;
- `pars::ProcessPars`: process parameters;
- `intorder::Int`: integration order;
- `elemtype::FiniteElement`: type of finite element.
"""
function calulate_stresses_3d_nonloc_depr(deformations::Array, displacements::Array, 
    elasticitymatrix::Matrix, beta_loc::Number, beta_nonloc::Number, neighbours::Array, 
    impactdist::Number, pars::ProcessPars, intorder::Int, elemtype::FiniteElement
)
    n_of_deformations_types = 6
    if size(deformations)[2] != n_of_deformations_types
        @error("Incorrect input deformations in calulate_stresses_3d_nonloc()")
        return nothing
    end

    # 1. Calculate local part of stresses
    locpart = calculateStresses(deformations, elasticitymatrix)
    locpart .*= beta_loc

    # 2. Calculate non-local part of stresses
    elements = pars.mesh.elements
    nodes = pars.mesh.nodes
    nonlocpart = zeros(length(nodes), n_of_deformations_types)
    averaging_nums = zeros(length(nodes))
    for elsrc_idx in eachindex(elements)
        elsrc_nodes = elements[elsrc_idx]
        nodes_per_elem_src = length(elsrc_nodes)
        x_source = [nodes[elsrc_nodes[i]][1] for i in 1:nodes_per_elem_src]
        y_source = [nodes[elsrc_nodes[i]][2] for i in 1:nodes_per_elem_src]
        z_source = [nodes[elsrc_nodes[i]][3] for i in 1:nodes_per_elem_src]

        for srcnodeidx in eachindex(elsrc_nodes)
            node_source = elsrc_nodes[srcnodeidx]

            loc_coords = getRSFromNode(srcnodeidx, elemtype)
            r_source = loc_coords[1]
            s_source = loc_coords[2]
            t_source = loc_coords[3]

            for elimp_idx in neighbours[elsrc_idx]
                elimp_nodes = elements[elimp_idx]
                nodes_per_elem_imp = length(elimp_nodes)
                x_impact = [nodes[elimp_nodes[i]][1] for i in 1:nodes_per_elem_imp]
                y_impact = [nodes[elimp_nodes[i]][2] for i in 1:nodes_per_elem_imp]
                z_impact = [nodes[elimp_nodes[i]][3] for i in 1:nodes_per_elem_imp]

                element_imp_displacements = []
                for node in elimp_nodes
                    push!(element_imp_displacements, displacements[3 * node - 2])
                    push!(element_imp_displacements, displacements[3 * node - 1])
                    push!(element_imp_displacements, displacements[3 * node])
                end

                f_integrand(r_imp, s_imp, t_imp) = nonloc_stresses_integrand(r_source,
                    s_source, t_source, r_imp, s_imp, t_imp, x_source, y_source, z_source,
                    x_impact, y_impact, z_impact, element_imp_displacements, impactdist,
                    elasticitymatrix, elemtype)

                defimpact = MultipleIntegral.gaussmethod_matrix_3d_nonsym(f_integrand, 
                    intorder)

                for typeidx in 1:n_of_deformations_types
                    nonlocpart[node_source, typeidx] += defimpact[typeidx]
                end
                averaging_nums[node_source] += 1
            end
        end
    end
    for nodeidx in 1:size(nonlocpart)[1]
        for typeidx in 1:size(nonlocpart)[2]
            if (averaging_nums[nodeidx] == 0)
                averaging_nums[nodeidx] = 1
            end
            nonlocpart[nodeidx, typeidx] /= averaging_nums[nodeidx]
        end
    end
    nonlocpart .*= beta_nonloc
    @show(nonlocpart)

    return locpart + nonlocpart
end  # calulate_stresses_nonloc

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

"""
    calculate_stresses_3d_nl(displ::Vector, pars::ProcessPars, beta_loc::Number, 
        beta_nl::Number, imp_dist::Number, elas_matr::Array, neighbours::Dict, 
        intorder::Int, elemtype::FiniteElement)

Calculate stresses according to non-local elasticity theory.

# Arguments
- `displ::Vector`: displacements vector;
- `pars::ProcessPars`: process parameters;
- `beta_loc::Number`: local part impact coefficient;
- `beta_nl::Number`: non-local part impact coefficient;
- `imp_dist::Number`: impact distance;
- `elas_matr::Array`: elasticity matrix;
- `neighbours::Dict`: global list of neighbours;
- `intorder::Int`: integration order;
- `elemtype::FiniteElement`: type of finite element.
"""
function calculate_stresses_3d_nl(displ::Vector, pars::ProcessPars, beta_loc::Number, 
    beta_nl::Number, imp_dist::Number, elas_matr::Array, neighbours::Dict, 
    intorder::Int, elemtype::FiniteElement
)
    n_nodes = length(displ)
    n_of_stresses_types = 6  # TODO: Compute this value
    stresses = zeros(Float64, n_nodes, n_of_stresses_types)
    els_cnts = zeros(n_nodes)

    for elemix in eachindex(pars.mesh.elements)
        element = pars.mesh.elements[elemix]
        nodes_per_elem = length(element)
        displ_source = []
        for node in element
            push!(displ_source, displ[3 * node - 2])
            push!(displ_source, displ[3 * node - 1])
            push!(displ_source, displ[3 * node])
        end

        for nodeix in eachindex(element)
            # Local part
            loc_coords = getRSFromNode(nodeix, elemtype)
            x_coords = [pars.mesh.nodes[element[i]][1] for i in 1:nodes_per_elem]
            y_coords = [pars.mesh.nodes[element[i]][2] for i in 1:nodes_per_elem]
            z_coords = [pars.mesh.nodes[element[i]][3] for i in 1:nodes_per_elem]

            grad_matr_source = gradMatr(loc_coords[1], loc_coords[2], loc_coords[3], 
                x_coords, y_coords, z_coords, elemtype)

            stresses_loc = beta_loc .* (elas_matr * grad_matr_source * displ_source)

            sourcept_glob = conv_loc_to_glob(loc_coords[1], loc_coords[2], loc_coords[3], 
                x_coords, y_coords, z_coords, elemtype)

            # Non-local part
            stresses_nl = zeros(length(stresses_loc))
            for elem_impact in neighbours[elemix]
                elem_nodes = pars.mesh.elements[elem_impact]
                x_coords_imp = [pars.mesh.nodes[elem_nodes[i]][1]
                    for i in 1:nodes_per_elem]
                y_coords_imp = [pars.mesh.nodes[elem_nodes[i]][2]
                    for i in 1:nodes_per_elem]
                z_coords_imp = [pars.mesh.nodes[elem_nodes[i]][3]
                    for i in 1:nodes_per_elem]

                displ_imp = []
                for node_imp in elem_nodes
                    push!(displ_imp, displ[3 * node_imp - 2])
                    push!(displ_imp, displ[3 * node_imp - 1])
                    push!(displ_imp, displ[3 * node_imp])
                end

                f_integrand(r_imp, s_imp, t_imp) = stresses_nl_integrand_3d(r_imp, s_imp, 
                    t_imp, x_coords_imp, y_coords_imp, z_coords_imp, elas_matr, imp_dist,
                    sourcept_glob, displ_imp, elemtype)

                stress = MultipleIntegral.gaussmethod_matrix_3d_nonsym(f_integrand, 
                    intorder)
                stresses_nl += stress
            end

            stresses_nl .*= beta_nl

            stresses_total = stresses_loc + stresses_nl

            node = element[nodeix]

            for direction in 1:n_of_stresses_types
                stresses[node, direction] += stresses_total[direction]
            end
            els_cnts[node] += 1
        end
    end

    for nodeidx in 1:size(stresses)[1]
        for typeidx in 1:size(stresses)[2]
            stresses[nodeidx, typeidx] /= els_cnts[nodeidx]
        end
    end

    return stresses
end

function stresses_nl_integrand_3d(r, s, t, x_coords, y_coords, z_coords, elas_matr::Array,
    imp_dist::Number, sourcept_glob::Tuple, elem_displ::Vector, elemtype::FiniteElement
)
    b_matr = gradMatr(r, s, t, x_coords, y_coords, z_coords, elemtype)

    # Calculating impact function
    impactpt_glob = conv_loc_to_glob(r, s, t, x_coords, y_coords, z_coords, elemtype)
    impactvec = vecfrompoints_3d(sourcept_glob, impactpt_glob)
    currdist = veclength_3d(impactvec)
    normfact = 1 / (8 * pi * imp_dist^3)
    impact = nonloc_gaussimpact(normfact, imp_dist, currdist)

    res_matr = impact .* (elas_matr * b_matr * elem_displ)

    return res_matr
end
