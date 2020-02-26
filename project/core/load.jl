# Loads and boundary conditions description

module Mod_Load

export bc

@enum bc begin
    fixedX  # Dislpacement fixed by X
    fixedY  # Dislpacement fixed by Y
    fixedXY  # Dislpacement fixed by X and Y
end

end  # Mod_Load