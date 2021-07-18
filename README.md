## Core

This branch describes the core of FEM application.

Basically FEM core is represented as Julia module named `CoreFEM`. It conains all FEM base methods that are required for appropriate calculations.

To make it work properly you must also import other FEM modules:
1. `BaseInterface` - to output calculation results;
2. `MeshFEM` - to use FEM mesh structure;
3. `MultipleIntegral` - to use integration methods.

Also you may want to use some stock finite element models for your calculations. For now it's only one model: `Quad4Pts`.
