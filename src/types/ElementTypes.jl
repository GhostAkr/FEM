module ElementTypes

export FiniteElement

export jacGlobToLoc, DetJs, gradMatr, displInterpMatr, nodesFromDirection, directionFromNodes, getRSFromNode

# TODO: Write macro to export whole enum at once
export FETypes, Quad4TypeID, Quad8TypeID

@enum FETypes begin
    Quad4TypeID
    Quad8TypeID
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

end