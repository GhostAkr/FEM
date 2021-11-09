# Process models

Model in current meaning is basically main expressions that comprehensively describe current finite element model. 
On this page you can learn how to use stock models and how to make your own.

## Guide

Each model is represented as Julia module, so you can just import your own one by native Julia way. Some stock models are already implemented and can be used
in your calculations. Now let's desrcibe how to make your own models.

To make your model you need to create several main methods in appropriate module:
1. `jacGlobToLoc(r, s, xCoords::Array{Float64}, yCoords::Array{Float64})` - function which returns Jacobi matrix for conversion between different coordinate systems;
2. `jacGlobToLocInv(r, s, xCoords::Array{Float64}, yCoords::Array{Float64})` - inversed Jacobi matrix;
3. `gradMatr(r, s, xCoords::Array{Float64}, yCoords::Array{Float64})` - gradient matrix of appropriate finite element model;
4. `displInterpMatr(r, s)` - displacements interpolation matrix of appropriate finite element model.

Note: it's important to use same names for appropriate functions above. And in fact it doesn't matter how these functions look inside.

You also can see examples in some stock models ([Stock models](@ref)).

## Stock models

### Quad4Pts

```@autodocs
Modules = [Quad4Pts]
Order = [:module, :type, :function]
```