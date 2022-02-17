# Examples

This page contains description of examples of FEM project which can be found in `examples`
folder. To run example you can use appropriate block of code in `Examples.jl`.

## SmallPlate

2D square plate, ``100 \times 100`` mm. Mesh in this case is pretty simple: 4 
evenly distributed elements.

![SmallPlate example](./images/smallplate.png)

### Tasks

- **Task**. Plate is fixed at one side and is stretched at the opposite one. Type of finite 
    element model: plain stress.

## Beam

2D beam, ``20 \times 100`` mm. Mesh consists of 2000 elements which are evenly distributed
over the beam.

![Beam example](./images/beam.png)

### Tasks

- **Task**. Beam is fixed at one side and is stretched at the opposite one. Type of finite 
    element model: plain stress.
