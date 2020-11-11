# TODO: Write macro to export whole enum at once
export meshType, Quad4Pts2D, Quad8Pts2D

"""
    meshType

Enum containing different mesh types:
- Quad4Pts2D;
- Quad8Pts2D;
"""
@enum meshType begin
    Quad4Pts2D
    Quad8Pts2D
end