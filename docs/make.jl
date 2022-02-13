# Documentation builder

# Set paths
include("../Set_Path.jl")

# List of all modules that should be included in documentation
using Documenter, MeshFEM, CoreFEM, Quad4Pts, MultipleIntegral

makedocs(sitename = "FEM Documentation")