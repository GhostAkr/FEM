using VectorFEM
using Iso8Pts3D
using Quad8Pts
using Quad4Pts
using MultipleIntegral
using LinearAlgebra

"""
    get_elem_neighbours!(neighbours::Array, elemnum::Int, distance::Number, 
        startpt::Tuple{Number, Number, Number}, parameters::ProcessPars)

Find all neighbours of given element. Neighbours are elements which are located not
more then `distance` from given one. This method calls itself recursively.

# Arguments
- `neighbours::Array`: numbers of neighbour elements;
- `elemnum::Int`: number of element for which neighbours should be found;
- `distance::Number`: distance which should be used to find neighbours;
- `startpt::Tuple{Number, Number, Number}: start point from which distance should be 
calculated;
- `parameters::ProcessPars`: simulation parameters.
"""
function get_elem_neighbours!(neighbours::Array, elemnum::Int, distance::Number, 
    startpt::Tuple{Number, Number, Number}, parameters::ProcessPars, 
    elem_type::FiniteElement
)
    nodes = parameters.mesh.elements[elemnum]
    for node in nodes
        # KLUDGE: due to historical restrictions one node can be represented as 2 
        # coordinates (x and y) instead of 3.
        coords_cnt = length(parameters.mesh.nodes[node])
        if coords_cnt == 2
            node_coords = (parameters.mesh.nodes[node][1], 
                parameters.mesh.nodes[node][2], 0)
        elseif coords_cnt == 3
            node_coords = parameters.mesh.nodes[node]
        else
            @error("Bad node in get_elem_neighbours!()")
            return nothing
        end

        # Find all element containing current node
        for elemidx in eachindex(parameters.mesh.elements)
            elem_nodes = parameters.mesh.elements[elemidx]
            if node in elem_nodes
                if elemidx in neighbours
                    continue
                end

                # Find element's center
                x = Vector{Float64}()
                y = Vector{Float64}()
                z = Vector{Float64}()
                for n in elem_nodes
                    append!(x, float(parameters.mesh.nodes[n][1]))
                    append!(y, float(parameters.mesh.nodes[n][2]))

                    if coords_cnt == 3
                        append!(z, float(parameters.mesh.nodes[n][3]))
                    else
                        append!(z, 0.)
                    end
                end
                glob_coords = conv_loc_to_glob(0., 0., 0., x, y, z, elem_type)

                # Check if element is inside impact area
                node_vec = vecfrompoints_3d(startpt, glob_coords)
                node_len = veclength_3d(node_vec)

                if node_len > distance
                    continue
                end

                append!(neighbours, elemidx)
                get_elem_neighbours!(neighbours, elemidx, distance, startpt, parameters,
                    elem_type)
            end
        end
    end
end

"""
	nonloc_gaussimpact(normfactor::Number, impactdistance::Number, distance::Number)

Computes non-local impact of point which is located at a distance `distance` from start 
one. Impact functions is represented as gauss function: ``a(r) = k exp(-r^2 / l^2)``, 
where ``k`` is normalization factor which can be obtained by solving normalization 
condition: ``\\int \\limits_{\\mathbb{R}^3} a(|x' - x|) \\text{d} V' = 1``, where
``r = |x' - x|``.

# Arguments
- `normfactor::Number`: normalization factor;
- `impactdistance::Number`: impact distance ``l``;
- `distance::Number`: distance from start point to current one.
"""
function nonloc_gaussimpact(normfactor::Number, impactdistance::Number, distance::Number)
    return normfactor * exp(-distance / impactdistance)
end

"""
	contribute_leftpart_nonloc!(pars::ProcessPars, targetmatr::Array, currmatr::Array, 
    source_elemnum::Int, impact_elemnum::Int, freedom_deg::Int)

Process a contribute of non-local stiffness matrix `currmatr` to global stiffness matrix
`targetmatr`.

# Arguments
- `pars::ProcessPars`: parameters of current process;
- `targetmatr::Array`: global stiffness matrix;
- `currmatr::Array`: non-local stiffness matrix;
- `source_elemnum::Int`: number of element which is under impact;
- `impact_elemnum::Int`: number of element which impacts;
- `freedom_deg::Int`: number of degrees of freedom.
"""
function contribute_leftpart_nonloc!(pars::ProcessPars, targetmatr::Array, currmatr::Array, 
	source_elemnum::Int, impact_elemnum::Int, freedom_deg::Int
)
    source_nodes = pars.mesh.elements[source_elemnum]
    impact_nodes = pars.mesh.elements[impact_elemnum]

    for i in eachindex(source_nodes)
        for j in eachindex(impact_nodes)
            for first_offset in 1:freedom_deg
                for second_offset in 1:freedom_deg
                    targetmatr[freedom_deg * source_nodes[i] - (first_offset - 1), 
                        freedom_deg * impact_nodes[j] - (second_offset - 1)] += 
                        currmatr[freedom_deg * i - (first_offset - 1), 
                        freedom_deg * j - (second_offset - 1)]
                end
            end
        end
    end
end

"""
	contribute_leftpart_by_connmatr_nonloc!(targetmatr::Array, 
    currmatr::Array, source_connmatr::Array, impact_connmatr::Array)

Process a contribute of non-local stiffness matrix `currmatr` to global stiffness matrix
`targetmatr` using given connectivity matrices.

# Arguments
- `targetmatr::Array`: global stiffness matrix;
- `currmatr::Array`: non-local stiffness matrix;
- `source_connmatr::Array`: connectivity matrix of source element;
- `impact_connmatr::Array`: connectivity matrix of impact element.
"""
function contribute_leftpart_by_connmatr_nonloc!(targetmatr::Array, currmatr::Array, 
    source_connmatr::Array, impact_connmatr::Array
)
    contribution = transpose(source_connmatr) * currmatr * impact_connmatr
    targetmatr .+= contribution
end

"""
    stiffnessintegrand_2d_nonloc(r_source, s_source, r_impact, s_impact,
        x_source::Array{Float64}, y_source::Array{Float64}, 
        x_impact::Array{Float64}, y_impact::Array{Float64},
        impactdist::Number, elasticitymatrix::AbstractArray, elemtype::FiniteElement)

Computes intgrand to calculate non-local stiffness matrix: ``F = A(r) B_e^T C B_{e'}``. 

# Arguments
- `r_source`: r-coordinate of node in source element;
- `s_source`: s-coordinate of node in source element;
- `r_impact`: r-coordinate of node in impact element;
- `s_impact`: s-coordinate of node in impact element;
- `x_source::Array{Float64}`: x-coordinates of nodes in source element;
- `y_source::Array{Float64}`: y-coordinates of nodes in source element;
- `x_impact::Array{Float64}`: x-coordinates of nodes in impact element;
- `y_impact::Array{Float64}`: y-coordinates of nodes in impact element;
- `impactdist::Number`: impact distance;
- `elasticitymatrix::AbstractArray`: elasticity matrix;
- `elemtype::FiniteElement`: type of finite element.
"""
function stiffnessintegrand_2d_nonloc(r_source, s_source, r_impact, s_impact, 
    x_source::Array{Float64}, y_source::Array{Float64}, 
    x_impact::Array{Float64}, y_impact::Array{Float64},
    impactdist::Number, elasticitymatrix::AbstractArray, elemtype::FiniteElement
)
    b_source = transpose(gradMatr(r_source, s_source, x_source, y_source, elemtype))
    b_impact = gradMatr(r_impact, s_impact, x_impact, y_impact, elemtype)
    jac_source = jacGlobToLoc(r_source, s_source, x_source, y_source, elemtype)
    jac_impact = jacGlobToLoc(r_impact, s_impact, x_impact, y_impact, elemtype)

    integrmatr = b_source * elasticitymatrix * b_impact * det(jac_source) * det(jac_impact)

    # Calculating impact function
    sourcept_glob = conv_loc_to_glob(r_source, s_source, x_source, y_source, elemtype)
    impactpt_glob = conv_loc_to_glob(r_impact, s_impact, x_impact, y_impact, elemtype)
    impactvec = vecfrompoints_2d(sourcept_glob, impactpt_glob)
    currdist = veclength_3d(impactvec)
    normfact = 1 / (2 * pi * impactdist^2)
    impact = nonloc_gaussimpact(normfact, impactdist, currdist)

    integrmatr .*= impact

    return integrmatr
end

"""
    stiffnessintegrand_3d_nonloc(r_source, s_source, t_source, 
        r_impact, s_impact, t_impact,
        x_source::Array{Float64}, y_source::Array{Float64}, z_source::Array{Float64}, 
        x_impact::Array{Float64}, y_impact::Array{Float64}, z_impact::Array{Float64},
        impactdist::Number, elasticitymatrix::AbstractArray, elemtype::FiniteElement)

Computes intgrand to calculate non-local stiffness matrix: ``F = A(r) B_e^T C B_{e'}``. 

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
- `impactdist::Number`: impact distance;
- `elasticitymatrix::AbstractArray`: elasticity matrix;
- `elemtype::FiniteElement`: type of finite element.
"""
function stiffnessintegrand_3d_nonloc(r_source, s_source, t_source, 
    r_impact, s_impact, t_impact, 
    x_source::Array{Float64}, y_source::Array{Float64}, z_source::Array{Float64}, 
    x_impact::Array{Float64}, y_impact::Array{Float64}, z_impact::Array{Float64},
    impactdist::Number, elasticitymatrix::AbstractArray, elemtype::FiniteElement
)
    # integr_terms_time = @elapsed begin
    # b_matrix_time = @elapsed begin
    b_src_time = @elapsed begin
    b_source = transpose(gradMatr(r_source, s_source, t_source, x_source, y_source, 
        z_source, elemtype))
    end
    b_imp_time = @elapsed begin
    b_impact = gradMatr(r_impact, s_impact, t_impact, x_impact, y_impact, z_impact, 
        elemtype)
    end
    # jac_time = @elapsed begin
    jac_source = jacGlobToLoc(r_source, s_source, t_source, x_source, y_source, z_source, 
        elemtype)
    jac_impact = jacGlobToLoc(r_impact, s_impact, t_impact, x_impact, y_impact, z_impact, 
        elemtype)
    # end

    # integr_calc_time = @elapsed begin
    integrmatr = b_source * elasticitymatrix * b_impact * det(jac_source) * det(jac_impact)
    # end

    # Calculating impact function
    # impact_calc_time = @elapsed begin
    sourcept_glob = conv_loc_to_glob(r_source, s_source, t_source, x_source, y_source, 
        z_source, elemtype)
    impactpt_glob = conv_loc_to_glob(r_impact, s_impact, t_impact, x_impact, y_impact, 
        z_impact, elemtype)
    impactvec = vecfrompoints_3d(sourcept_glob, impactpt_glob)
    currdist = veclength_3d(impactvec)
    normfact = 1 / (8 * pi * impactdist^3)
    impact = nonloc_gaussimpact(normfact, impactdist, currdist)
    # end

    # impact_contr_time = @elapsed begin
    integrmatr .*= impact
    # end

    # println("Calculating terms for integrand: ", integr_terms_time)
    # println("Calculating integrand: ", integr_calc_time)
    # println("Calculating impact function: ", impact_calc_time)
    # println("Impact function contribution: ", impact_contr_time)
    # println("B-matrix calculation: ", b_matrix_time)
    # println("Jac calculation: ", jac_time)
    # println("Source B time: ", b_src_time)
    # println("Impact B time: ", b_imp_time)
    # println()

    return integrmatr
end

"""
    stiffnessmatr_2d_nonloc(elasmatr::Array, parameters::ProcessPars, 
        source_elemnum::Int, impact_elemnum::Int, impactdist::Number, intorder::Int, 
        elemtype::FiniteElement)

Calcuulate non-local stiffness matrix for element number `source_elemnum` with 
respect to impact of element number `impact_elemnum`.

# Arguments
- `elasmatr::Array`: elasticity matrix;
- `parameters::ProcessPars`: parameters of current process;
- `source_elemnum::Int`: number of element which is under impact;
- `impact_elemnum::Int`: number of element which impacts;
- `impactdist::Number`: impact distance;
- `intorder::Int`: integration order;
- `elemtype::FiniteElement`: type of finite element.
"""
function stiffnessmatr_2d_nonloc(elasmatr::Array, parameters::ProcessPars, 
	source_elemnum::Int, impact_elemnum::Int, impactdist::Number, intorder::Int, 
	elemtype::FiniteElement
)
    nodes_source_cnt = length(parameters.mesh.elements[source_elemnum])
    nodes_impact_cnt = length(parameters.mesh.elements[impact_elemnum])

    x_source = [parameters.mesh.nodes[parameters.mesh.elements[source_elemnum][i]][1] 
        for i in 1:nodes_source_cnt]
    y_source = [parameters.mesh.nodes[parameters.mesh.elements[source_elemnum][i]][2] 
        for i in 1:nodes_source_cnt]

    x_impact = [parameters.mesh.nodes[parameters.mesh.elements[impact_elemnum][i]][1] 
        for i in 1:nodes_impact_cnt]
    y_impact = [parameters.mesh.nodes[parameters.mesh.elements[impact_elemnum][i]][2] 
        for i in 1:nodes_impact_cnt]

    f_integrand(r_loc, s_loc, r_imp, s_imp) = stiffnessintegrand_2d_nonloc(r_loc, s_loc, 
        r_imp, s_imp, x_source, y_source, x_impact, y_impact, impactdist, elasmatr, elemtype
    )

    nonloc_matr = MultipleIntegral.gaussmethod_matrix_2d_nonloc(f_integrand, intorder)
    return nonloc_matr
end

"""
    stiffnessmatr_3d_nonloc(elasmatr::Array, parameters::ProcessPars, 
        source_elemnum::Int, impact_elemnum::Int, impactdist::Number, intorder::Int, 
        elemtype::FiniteElement)

Calcuulate non-local stiffness matrix for element number `source_elemnum` with 
respect to impact of element number `impact_elemnum`.

# Arguments
- `elasmatr::Array`: elasticity matrix;
- `parameters::ProcessPars`: parameters of current process;
- `source_elemnum::Int`: number of element which is under impact;
- `impact_elemnum::Int`: number of element which impacts;
- `impactdist::Number`: impact distance;
- `intorder::Int`: integration order;
- `elemtype::FiniteElement`: type of finite element.
"""
function stiffnessmatr_3d_nonloc(elasmatr::Array, parameters::ProcessPars, 
	source_elemnum::Int, impact_elemnum::Int, impactdist::Number, intorder::Int, 
	elemtype::FiniteElement
)
    # @time begin
    nodes_source_cnt = length(parameters.mesh.elements[source_elemnum])
    nodes_impact_cnt = length(parameters.mesh.elements[impact_elemnum])

    x_source = [parameters.mesh.nodes[parameters.mesh.elements[source_elemnum][i]][1] 
        for i in 1:nodes_source_cnt]
    y_source = [parameters.mesh.nodes[parameters.mesh.elements[source_elemnum][i]][2] 
        for i in 1:nodes_source_cnt]
    z_source = [parameters.mesh.nodes[parameters.mesh.elements[source_elemnum][i]][3] 
        for i in 1:nodes_source_cnt]

    x_impact = [parameters.mesh.nodes[parameters.mesh.elements[impact_elemnum][i]][1] 
        for i in 1:nodes_impact_cnt]
    y_impact = [parameters.mesh.nodes[parameters.mesh.elements[impact_elemnum][i]][2] 
        for i in 1:nodes_impact_cnt]
    z_impact = [parameters.mesh.nodes[parameters.mesh.elements[impact_elemnum][i]][3] 
        for i in 1:nodes_impact_cnt]
    # end

    # @time begin
    f_integrand(r_loc, s_loc, t_loc, r_imp, s_imp, t_imp) = stiffnessintegrand_3d_nonloc(
        r_loc, s_loc, t_loc, r_imp, s_imp, t_imp, x_source, y_source, z_source, x_impact, 
        y_impact, z_impact, impactdist, elasmatr, elemtype
    )
    # end

    # @time begin
    nonloc_matr = MultipleIntegral.gaussmethod_matrix_3d_nonloc(f_integrand, intorder)
    # end
    return nonloc_matr
end
