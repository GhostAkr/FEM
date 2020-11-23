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

"""
    loadDirection

Enum containing different load local directions for isoparametric 
finite element. Considering native system ``(r, s)`` where ``r`` is directed to the right
and ``s`` is directed to the top.
- `right`: ``+r`` direction;
- `left`: ``-r`` direction;
- `top`: ``+s`` direction;
- `bottom`: ``-s`` direction.
"""
@enum loadDirection begin
    top = 1
    left = 2
    bottom = 3
    right = 4
    backwards = 5  # To us
    towards = 6  # From us
end
