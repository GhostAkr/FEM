# Material properties

module Mod_Mat_Props

    # Possible material properties
    @enum materialProperty begin
        poisC  # Poisson's coefficient
        youngMod  # Young's modulus
    end

    # Export block
    export materialProperty

end  # Mod_Mat_Props