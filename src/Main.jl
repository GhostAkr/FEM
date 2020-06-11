# To provide correct modules loading it's important to execeute Set_Path.jl first.

using CoreFEM
using multipleIntegral
using MeshFEM
using BaseInterface

Res = CoreFEM.fem2D("examples/Plate_Mesh.dat", "examples/Input")
