# Material properties

module Material

export materialProperty

# Possible material properties
@enum materialProperty begin
    poisC  # Poisson's coefficient
    youngMod  # Young's modulus
end

end  # Material