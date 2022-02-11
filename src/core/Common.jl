"""
    get_connmatr(pars::processPars, elemnum::Int, freedom_deg::Int)

Computes connectivity matrix for nodes in given element.

# Arguments
- `pars::processPars`: parameters of current process;
- `elemnum::Int`: number of finite element in mesh;
- `freedom_deg::Int`: number of degrees of freedom.
"""
function get_connmatr(pars::processPars, elemnum::Int, freedom_deg::Int)
    elem_nodes = pars.mesh.elements[elemnum]
    nodes_cnt = length(pars.mesh.nodes)

    rows_cnt = length(elem_nodes) * freedom_deg
    cols_cnt = nodes_cnt * freedom_deg
    connmatr = zeros(rows_cnt, cols_cnt)

    for node in eachindex(elem_nodes)
        for dof in 1:freedom_deg
            connmatr[node * freedom_deg - (dof - 1), elem_nodes[node] * freedom_deg - 
                (dof - 1)] = 1
        end
    end

    return connmatr
end