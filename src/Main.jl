# Entry point of application. Run this file in Julia REPL to load all necessary variables
# to current session.

# Adding paths for current session to properly load all local Julia modules.
include("../Set_Path.jl")

# Modules which are used by this application
using CoreFEM           # Contains all computational routine
using ElementTypes      # Contains information about finite elements type system
