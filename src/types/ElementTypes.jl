module ElementTypes

export FiniteElement

export jacGlobToLoc, DetJs, gradMatr, displInterpMatr, nodesFromDirection, getRSFromNode

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

DetJs(r, s, xCoords::Array{Float64}, yCoords::Array{Float64}, elType::FiniteElement) = nothing

gradMatr(r, s, xCoords::Array{Float64}, yCoords::Array{Float64}, elType::FiniteElement) = nothing

displInterpMatr(r, s) = nothing

nodesFromDirection(direction::Int, elType::FiniteElement) = nothing

getRSFromNode(nodeIndex::Int) = nothing

end