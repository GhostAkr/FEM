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

## Beam3D. Small task
3D beam, ``10 \times 10 \times 100`` mm. Mesh consists of 40 elements. On the end of beam 
mesh refinement is presented.

![Beam3D/Small example](./images/beam3D_small.png)

### Tasks

- **TaskStretch**. Beam is fixed at one side. On the opposite side pressure is applied to
    the end of beam.
- **TaskStretchNL**. The same as usual **TaskStretch** but with non-local parameters:
    - local impact: ``0.8``;
    - non-local impact: ``0.2``;
    - impact distance: ``10``.


## Beam3D. Big task

3D beam, ``10 \times 10 \times 100`` mm. Mesh consists of 2500 elements which are evenly
distributed over the beam.

![Beam3D/Big example](./images/beam3D_big.png)

### Tasks

- **TaskBind**. Beam is fixed at one side. On the oppsite side pressure is applied to the 
    top of beam.

## Beam3D. Analogue 2D

3D beam, ``10 \times 10 \times 100`` mm. Mesh consists only of 14 elements. Main feature of
this example is that mesh has only one element along Z axis.

![Beam3D/Analogue 2D](./images/beam3D_analogue2d.png)

### Tasks

- **TaskStretch**. Beam is fixed at one side. On the opposite side pressure is applied to
    the end of beam.
- **TaskStretchNL**. The same as usual **TaskStretch** but with non-local parameters:
    - local impact: ``0.8``;
    - non-local impact: ``0.2``;
    - impact distance: ``10``.

## Beam3D. Mini task

3D beam, ``10 \times 10 \times 100`` mm. Mesh consists only of 16 elements. Initially this
example was created to investigate non-local model behaviour. For example we can easily
check how neighbours interaction works for each element.

![Beam3D/Mini example](./images/beam3D_mini.png)

### Tasks

- **TaskStretch**. Beam is fixed at one side. On the opposite side pressure is applied to
    the end of beam.
- **TaskStretchNL**. The same as usual **TaskStretch** but with non-local parameters:
    - local impact: ``0.8``;
    - non-local impact: ``0.2``;
    - impact distance: ``10``.
