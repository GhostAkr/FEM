using VectorFEM

"""
    get_elem_neighbours!(neighbours::Array, elemnum::Int, distance::Number, 
        startpt::Tuple{Number, Number, Number}, parameters::processPars)

Find all neighbours of given element. Neighbours are elements which are located not
more then `distance` from given one. This method calls itself recursively.

# Arguments
- `neighbours::Array`: numbers of neighbour elements;
- `elemnum::Int`: number of element for which neighbours should be found;
- `distance::Number`: distance which should be used to find neighbours;
- `startpt::Tuple{Number, Number, Number}: start point from which distance should be 
calculated;
- `parameters::processPars`::
"""
function get_elem_neighbours!(neighbours::Array, elemnum::Int, distance::Number, 
    startpt::Tuple{Number, Number, Number}, parameters::processPars
)
    nodes = parameters.mesh.elements[elemnum]
    for node in nodes
        # KLUDGE: due to historical restrictions one node can be represented as 2 
        # coordinates (x and y) instead of 3.
        nodes_cnt = length(parameters.mesh.nodes[node + 1])
        if nodes_cnt == 2
            node_coords = (parameters.mesh.nodes[node + 1][1], 
                parameters.mesh.nodes[node + 1][2], 0)
        elseif nodes_cnt == 3
            node_coords = parameters.mesh.nodes[node + 1]
        else
            @error("Bad node in get_elem_neighbours!()")
            return nothing
        end

        node_vec = vecfrompoints_3d(startpt, node_coords)
        node_len = veclength_3d(node_vec)

        if node_len > distance
            continue
        end

        # Find all element containing current node
        for elemidx in eachindex(parameters.mesh.elements)
            elem = parameters.mesh.elements[elemidx]
            if node in elem
                if elemidx in neighbours
                    continue
                end
                append!(neighbours, elemidx)
                get_elem_neighbours!(neighbours, elemidx, distance, startpt, parameters)
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
    return normfactor * exp(-distance^2 / impactdistance^2)
end

"""
	contribute_leftpart_nonloc!(pars::processPars, targetmatr::Array, currmatr::Array, 
    source_elemnum::Int, impact_elemnum::Int, freedom_deg::Int)

Process a contribute of non-local stiffness matrix `currmatr` to global stiffness matrix
`targetmatr`.

# Arguments
- `pars::processPars`: parameters of current process;
- `targetmatr::Array`: global stiffness matrix;
- `currmatr::Array`: non-local stiffness matrix;
- `source_elemnum::Int`: number of element which is under impact;
- `impact_elemnum::Int`: number of element which impacts;
- `freedom_deg::Int`: number of degrees of freedom.
"""
function contribute_leftpart_nonloc!(pars::processPars, targetmatr::Array, currmatr::Array, 
	source_elemnum::Int, impact_elemnum::Int, freedom_deg::Int
)
    source_nodes = pars.mesh.elements[source_elemnum]
    impact_nodes = pars.mesh.elements[impact_elemnum]

    for i in eachindex(source_nodes)
        for j in eachindex(impact_nodes)
            for first_offset in 1:freedom_deg
                for second_offset in 1:freedom_deg
                    targetMatrix[freedom_deg * source_nodes[i] - (first_offset - 1), 
                        freedom_deg * impact_nodes[j] - (second_offset - 1)] += 
                        currentElementMatrix[freedom_deg * i - (first_offset - 1), 
                        freedom_deg * j - (second_offset - 1)]
                end
            end
        end
    end
end
