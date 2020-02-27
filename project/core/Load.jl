# Loads and boundary conditions description

module Load

export bc

@enum bc begin
    fixedX  # Dislpacement fixed by X
    fixedY  # Dislpacement fixed by Y
    fixedXY  # Dislpacement fixed by X and Y
end

end  # Load