# Main file. Used as application start point

# Main module
module Mod_Main

    # Include block
    include("include.jl")
    include("core/fem.jl")

    # Import block
    import Base
    import Mod_FEM

    Mod_FEM.fem2D()
end  # Mod_Main
