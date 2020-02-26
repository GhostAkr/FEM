# Material properties

module Mod_Mat_Props

export materialProperty

# Possible material properties
@enum materialProperty begin
    poisC  # Poisson's coefficient
    youngMod  # Young's modulus
end

end  # Mod_Mat_Props