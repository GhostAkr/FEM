# Scripts for setting modules path

# Engine path
push!(LOAD_PATH, "src/mesh")
push!(LOAD_PATH, "src/core")
push!(LOAD_PATH, "src/core/math")
push!(LOAD_PATH, "src/core/interpolation")
push!(LOAD_PATH, "src/core/integration")
push!(LOAD_PATH, "src/interface")
push!(LOAD_PATH, "src")

# Documentation path
push!(LOAD_PATH, "docs/src/")

# Types path
push!(LOAD_PATH, "src/types")
