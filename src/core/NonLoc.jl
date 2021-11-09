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