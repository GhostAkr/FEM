# To provide correct modules loading it's important to execeute Set_Path.jl first.

using CoreFEM
using multipleIntegral
using MeshFEM
using BaseInterface

Res = CoreFEM.fem2D("examples/QuarterPlate/Quarter_Mesh.dat", "examples/QuarterPlate/Quarter_Data")
# Res = CoreFEM.fem2D("examples/QuarterPlate/Small_Mesh.dat", "examples/QuarterPlate/Small_Input")
