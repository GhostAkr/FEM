# Loads and boundary conditions description

module Mod_Load
    @enum bc begin
        fixedX  # Dislpacement fixed by X
        fixedY  # Dislpacement fixed by Y
        fixedXY  # Dislpacement fixed by X and Y
    end

    # Export block
    export bc
end  # Mod_Load