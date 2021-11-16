module ElementTypes

export FiniteElement

export jacGlobToLoc, DetJs, gradMatr, displInterpMatr, nodesFromDirection
export directionFromNodes, getRSFromNode, conv_loc_to_glob

# TODO: Write macro to export whole enum at once
export FETypes, Quad4TypeID, Quad8TypeID, Iso8Pts3DTypeID

@enum FETypes begin
    Quad4TypeID
    Quad8TypeID
    Iso8Pts3DTypeID
    # Don't forget to export new type if needed
end

abstract type FiniteElement end

jacGlobToLoc(r, s, xCoords::Array{Float64}, yCoords::Array{Float64}, elType::FiniteElement) = nothing

DetJs(r, s, xCoords::Array{Float64}, yCoords::Array{Float64}, load_direction::Int, elType::FiniteElement) = nothing

gradMatr(r, s, xCoords::Array{Float64}, yCoords::Array{Float64}, elType::FiniteElement) = nothing

displInterpMatr(r, s) = nothing

nodesFromDirection(direction::Int, elType::FiniteElement) = nothing
directionFromNodes(nodes::Array, elType:: FiniteElement) = nothing

getRSFromNode(nodeIndex::Int) = nothing

conv_loc_to_glob(r, s, x_coords::Array{Float64}, y_coords::Array{Float64}) = nothing

end