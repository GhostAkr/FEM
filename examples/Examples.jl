# File contains commands which correspond to preloaded examples. To run example you need to
# load all necessary modules by launching src/Main.jl file in Julia REPL. That will create
# all necessary variables in current Julia session. 

# 2D models

# SmallPlate example
begin
    CoreFEM.fem2D("examples/SmallPlate/Mesh.med", "examples/SmallPlate/Task.json", 
        Quad4TypeID)
end

# Beam example
begin
    CoreFEM.fem2D("examples/Beam/Mesh.med", "examples/Beam/Task.json", Quad4TypeID)
end

# 3D models

# Beam3D/SmallTask
begin
    CoreFEM.elasmech_3d("examples/Beam3D/SmallTask/Mesh.med",
        "examples/Beam3D/SmallTask/TaskStretch.json", Iso8Pts3DTypeID)
end

# Beam3D/SmallTask (non-local model)
begin
    CoreFEM.elasmech_3d_nonloc("examples/Beam3D/SmallTask/Mesh.med",
        "examples/Beam3D/SmallTask/TaskStretchNL.json", Iso8Pts3DTypeID)
end

# Beam3D/BigTask
begin
    CoreFEM.elasmech_3d("examples/Beam3D/BigTask/Mesh.med", 
        "examples/Beam3D/BigTask/TaskBind.json", Iso8Pts3DTypeID)
end

# Beam3D/Analogue2D
begin
    CoreFEM.elasmech_3d("examples/Beam3D/Analogue2D/Mesh.med",
        "examples/Beam3D/Analogue2D/TaskStretch.json", Iso8Pts3DTypeID)
end

# Beam3D/Analogue2D (non-local model)
begin
    CoreFEM.elasmech_3d_nonloc("examples/Beam3D/Analogue2D/Mesh.med", 
        "examples/Beam3D/Analogue2D/TaskStretchNL.json", Iso8Pts3DTypeID)
end

# Beam3D/MiniTask
begin
    CoreFEM.elasmech_3d("examples/Beam3D/MiniTask/Mesh.med", 
        "examples/Beam3D/MiniTask/TaskStretch.json", Iso8Pts3DTypeID)
end

# Beam3D/MiniTask (non-local model)
begin
    CoreFEM.elasmech_3d_nonloc("examples/Beam3D/MiniTask/Mesh.med", 
        "examples/Beam3D/MiniTask/TaskStretchNL.json", Iso8Pts3DTypeID)
end

# Bolt/MidTask
begin
    CoreFEM.elasmech_3d("examples/Bolt/MidTask/Mesh.med", 
        "examples/Bolt/MidTask/TaskStretch.json", Iso8Pts3DTypeID)
end

# Bolt/MidTask (non-local model)
begin
    CoreFEM.elasmech_3d_nonloc("examples/Bolt/MidTask/Mesh.med", 
        "examples/Bolt/MidTask/TaskStretchNL.json", Iso8Pts3DTypeID)
end
