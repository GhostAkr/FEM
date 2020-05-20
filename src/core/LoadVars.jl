# Base variables and types of load part

export bc

"""
    bc

Enum containing different boundary condition cases:
- `fixedX` - dislpacement fixed by X;
- `fixedY` - dislpacement fixed by Y;
- `fixedXY` - dislpacement fixed by X and Y.
"""
@enum bc begin
    fixedX
    fixedY
    fixedXY
end
