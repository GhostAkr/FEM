"""
    CoreFEM
Module describing FEM core.
"""
module CoreFEM

include("Input.jl")

using MeshFEM
using DelimitedFiles
using IterativeSolvers
using ElementTypes
using TestFEM
using BenchmarkTools
using SparseArrays
using Base.Threads

using Quad4Pts
using Quad8Pts
using Iso8Pts3D

include("MaterialVars.jl")
include("Load.jl")
include("StiffnessMatrix.jl")
include("Deformations.jl")
include("Constants.jl")
include("Stresses.jl")
include("NonLoc.jl")
include("Common.jl")

# Interface functions
include("interface/ExportResult.jl")
include("interface/ImportParams.jl")
include("interface/ImportJSON.jl")

export fem2D

"""
    assembly_left_part!(pars::ProcessPars, targetMatrix::Array, currentElementMatrix::Array, elementNum::Number, freedom_deg::Int)

Assemble left part of linear system of equations. This method applies given local stiffness matrix to global ensemble.

# Arguments
- `pars::ProcessPars`: parameters of current model;
- `targetMatrix::Array`: global stiffness matrix that should be updated;
- `currentElementMatrix::Array`: local stiffness matrix of given element that should be applid to global ensemble;
- `elementNum::Number`: number of given element in current mesh;
- `freedom_deg::Int`: degree of freedom.
"""
function assembly_left_part!(pars::ProcessPars, targetMatrix::Array, currentElementMatrix::Array, elementNum::Number, freedom_deg::Int)
    elementNodes = pars.mesh.elements[elementNum]
    for i in eachindex(elementNodes)
        for j in eachindex(elementNodes)
            for first_offset in 1:freedom_deg
                for second_offset in 1:freedom_deg
                    targetMatrix[freedom_deg * elementNodes[i] - (first_offset - 1), freedom_deg * elementNodes[j] - (second_offset - 1)] += 
                                currentElementMatrix[freedom_deg * i - (first_offset - 1), freedom_deg * j - (second_offset - 1)]
                end
            end
        end
    end
end  # assemblyFEM2D

function assemblyFEM3D(pars::ProcessPars, targetMatrix::Array, currentElementMatrix::Array, elementNum::Number)
    elementNodes = pars.mesh.elements[elementNum]
    for i in eachindex(elementNodes)
        for j in eachindex(elementNodes)
            # TODO: Fix current matrix indices
            targetMatrix[3 * elementNodes[i] - 2, 3 * elementNodes[j] - 2] += currentElementMatrix[3 * i - 2, 3 * j - 2]
            targetMatrix[3 * elementNodes[i] - 2, 3 * elementNodes[j] - 1] += currentElementMatrix[3 * i - 2, 3 * j - 1]
            targetMatrix[3 * elementNodes[i] - 2, 3 * elementNodes[j]] += currentElementMatrix[3 * i - 2, 3 * j]
            targetMatrix[3 * elementNodes[i] - 1, 3 * elementNodes[j] - 2] += currentElementMatrix[3 * i - 1, 3 * j - 2]
            targetMatrix[3 * elementNodes[i] - 1, 3 * elementNodes[j] - 1] += currentElementMatrix[3 * i - 1, 3 * j - 1]
            targetMatrix[3 * elementNodes[i] - 1, 3 * elementNodes[j]] += currentElementMatrix[3 * i - 1, 3 * j]
            targetMatrix[3 * elementNodes[i], 3 * elementNodes[j] - 2] += currentElementMatrix[3 * i, 3 * j - 2]
            targetMatrix[3 * elementNodes[i], 3 * elementNodes[j] - 1] += currentElementMatrix[3 * i, 3 * j - 1]
            targetMatrix[3 * elementNodes[i], 3 * elementNodes[j]] += currentElementMatrix[3 * i, 3 * j]
        end
    end
end  # assemblyFEM3D

"""
    applyConstraints(pars::ProcessPars, loads::Array, globalK::Array)

Applies constraints from model parameters to given ensemble.

# Arguments
- `pars::ProcessPars`: parameters of current model;
- `loads::Array`: global loads vector (right part of equations system);
- `globalK::Array`: global stiffness matrix (left part of equations system).
"""
function applyConstraints(pars::ProcessPars, loads::Array, globalK::Array)
    for (node, bc) in pars.bc
        if bc == fixedX
            applyFixedX(node, loads, globalK)
        elseif bc == fixedY
            applyFixedY(node, loads, globalK)
        elseif bc == fixedXY
            applyFixedXY(node, loads, globalK)
        else
            println("Unhandled boundary condition")
        end
    end
end  # applyConstraints

function applyConstraints3D(pars::ProcessPars, loads::Array, globalK::Array)
    for (node, bc) in pars.bc
        if bc == fixedX
            applyFixedX3D(node, loads, globalK)
        elseif bc == fixedY
            applyFixedY3D(node, loads, globalK)
        elseif bc == fixedXY
            applyFixedX3D(node, loads, globalK)
            applyFixedY3D(node, loads, globalK)
        elseif bc == fixedZ
            applyFixedZ3D(node, loads, globalK)
        elseif bc == fixedXZ
            applyFixedX3D(node, loads, globalK)
            applyFixedZ3D(node, loads, globalK)
        elseif bc == fixedYZ
            applyFixedY3D(node, loads, globalK)
            applyFixedZ3D(node, loads, globalK)
        elseif bc == fixedXYZ
            applyFixedX3D(node, loads, globalK)
            applyFixedY3D(node, loads, globalK)
            applyFixedZ3D(node, loads, globalK)
        else
            println("Unhandled boundary condition")
        end
    end
end  # applyConstraints

"""
    solve(globalK::Array, loadVector::Array)

Solve given equation system.

Since Julia provides native workaround with linear algebra elements the best solution in most cases
is to use standart syntax to solve linear equations system.
According to official documentation Julia will choose the best solving method by itself.
If it's not, there should be a way to control it.

# Arguments
- `globalK::Array`: global stiffness matrix (left part of equations system);
- `loadVector::Array`: global loads vector (right part of equations system).
"""
function solve(globalK::Array, loadVector::Array)
    return globalK \ loadVector
end

function defineElemType(elemTypeID::FETypes)
    resElement = nothing
    if elemTypeID === Quad4TypeID
        resElement = Quad4Type("Quad4Type")
    elseif elemTypeID === Quad8TypeID
        resElement = Quad8Type("Quad8Type")
    elseif elemTypeID === Iso8Pts3DTypeID
        resElement = Iso8Pts3DType("Iso8Pts3DType")
    else
        println("Unknown finite element type")
    end
    return resElement
end

function typeMeshFromElement(elemTypeID::FETypes)
    resMeshType = nothing
    if elemTypeID === Quad4TypeID
        resMeshType = Quad4Pts2D
    elseif elemTypeID === Quad8TypeID
        resMeshType = Quad8Pts2D
    elseif elemTypeID === Iso8Pts3DTypeID
        resMeshType = Iso8Pts3DMeshType
    else
        println("Unknown element type while converting to mesh type")
    end
    return resMeshType
end

"""
    fem2D()

Start calculation with given model.

# Arguments
- `meshPath::String`: Path to given mesh;
- `dataPath::String`: Path to given initial data.
"""
function fem2D(meshPath::String, dataPath::String, elemTypeID::FETypes)
    freedom_deg = 2
    # Getting element type
    elementType = defineElemType(elemTypeID)
    if (elementType === nothing)
        println("Element type passed to fem2D() is unknown")
        return
    end

    # Getting mesh type
    meshType = typeMeshFromElement(elemTypeID)

    # Model type
    model_type = default_modeltype

    parameters = ProcessPars(testMaterialProperties(), testBC(), testLoad(), 
        generateTestMesh2D(2), ModelPars(model_type))

    # Reading mesh
    split_path = splitext(meshPath)
    if (split_path[2] == ".dat")
        parameters.mesh = readMeshFromSalomeDAT(meshPath, meshType)
    elseif (split_path[2] == ".med")
        parameters.mesh = read_mesh_from_med(meshPath, meshType)
    else
        @error("Given mesh has unknown format")
        return nothing
    end

    # Reading parameters
    read_params_JSON!(dataPath, parameters)

    intOrder = 3
    nu = parameters.materialProperties[poisC]
    E = parameters.materialProperties[youngMod]
    C = elasticityMatrix(E, nu, parameters.model.type)
    ensembleMatrix = zeros(Float64, 2 * size(parameters.mesh.nodes)[1], 2 * size(parameters.mesh.nodes)[1])
    for elementNum in eachindex(parameters.mesh.elements)
        K = stiffnessMatrix(C, parameters, elementNum, intOrder, elementType)
        assembly_left_part!(parameters, ensembleMatrix, K, elementNum, freedom_deg)
    end
    loadVector = assembly_loads!(parameters, intOrder, elementType, freedom_deg)
    applyConstraints(parameters, loadVector, ensembleMatrix)

    mkpath("equation")
    # Writing left part to file
    open("equation/K", "w") do file
        writedlm(file, ensembleMatrix)
    end
    # Writing right part to file
    open("equation/F", "w") do file
        writedlm(file, loadVector)
    end
    println("Solving...")
    result = solve(ensembleMatrix, loadVector)

    if TestFEM.verify_example(meshPath, dataPath, result)
        @info "Result is correct"
    else
        @info "Result is INcorrect"
    end

    deformations = calculateDeformations(result, parameters, intOrder, elementType)
    stresses = calculateStresses(deformations, C)
    vonMises = calculateVonMises(stresses)
    # Writing result to file
    open("equation/result", "w") do file
        writedlm(file, result)
    end
    # Exporting results to VTK file
    exportToVTK(result, deformations, stresses, vonMises, parameters, meshType)
    return result
end  # fem2D

"""
    elasmech_3d(mesh_path::String, data_path::String, elem_type_id::FETypes)

Start calculation of 3D elastic mechanical problem.

# Arguments
- `mesh_path::String`: Path to mesh;
- `data_path::String`: Path to initial data;
- `elem_type_id::FETypes`: ID of finite element type.
"""
function elasmech_3d(mesh_path::String, data_path::String, elem_type_id::FETypes)
    freedom_deg = 3

    # 1. Getting element type
    element_type = defineElemType(elem_type_id)
    if (element_type === nothing)
        @error("Element type passed to fem2D() is unknown")
        return
    end

    # 2. Getting mesh type
    mesh_type = typeMeshFromElement(elem_type_id)

    # Model type
    model_type = default_modeltype

    parameters = ProcessPars(testMaterialProperties(), testBC3D(), testLoad3D(), 
        generateTestMesh3D(), ModelPars(model_type))

    # 3. Reading mesh
    parameters.mesh = read_mesh_from_med(mesh_path, mesh_type)

    # 4. Reading problem data
    read_params_JSON!(data_path, parameters)

    # 5. Integration order
    int_order = 2

    # 6. Problem constants
    nu = parameters.materialProperties[poisC]
    E = parameters.materialProperties[youngMod]
    C = elasticityMatrix(E, nu, parameters.model.type)

    println("Number of threads: ", nthreads())

    # 7. Stiffness matrix (left part of final equation)
    @info("Assembling left part...")
    @time begin
    ensemble_matrix = zeros(3 * size(parameters.mesh.nodes)[1], 3 * size(parameters.mesh.nodes)[1])
    @threads for element_num in eachindex(parameters.mesh.elements)
        k = stiffnessMatrix3D(C, parameters, element_num, int_order, element_type)
        assembly_left_part!(parameters, ensemble_matrix, k, element_num, freedom_deg)
    end
    end  # @time

    # 8. Load vector (right part og final equation)
    @info("Assembling right part...")
    @time begin
    load_vector = assembly_loads!(parameters, int_order, element_type, freedom_deg)
    end  # @time

    # 9. Applying constraints to equation
    @info("Applying constraints...")
    @time begin
    applyConstraints3D(parameters, load_vector, ensemble_matrix)
    end

    # 10. Creating folder for the equation entities
    #mkpath("equation")

    # 11. Writing left part to file
    #open("equation/K", "w") do file
    #    writedlm(file, ensemble_matrix)
    #end

    # 12. Writing right part to file
    #open("equation/F", "w") do file
    #    writedlm(file, load_vector)
    #end

    # 13. Solving equation
    @info("Solving...")
    @time begin
    result = solve(ensemble_matrix, load_vector)
    end  # @time

    # 14. Verifying result
    #if TestFEM.verify_example(mesh_path, data_path, result)
    #    @info "Result is correct"
    #else
    #    @info "Result is INcorrect"
    #end

    # 15. Writing result to file
    #open("equation/result", "w") do file
    #    writedlm(file, result)
    #end

    # 16. Calculating deformations
    deformations = calculate_deformations_3d(result, parameters, element_type)

    # 16. Exporting result to VTK
    exportToVTK(result, deformations, nothing, nothing, parameters, mesh_type)

    return result
end  # fem3D

"""
	elasmech_2d_nonloc(mesh_path::String, data_path::String, elem_type_id::FETypes)

Start calculation of 2D elastic non-local mechanical problem.

# Arguments
- `mesh_path::String`: path to mesh;
- `data_path::String`: path to initial data;
- `elem_type_id::FETypes`: ID of finite element type.
"""
function elasmech_2d_nonloc(mesh_path::String, data_path::String, elem_type_id::FETypes
)
    freedom_deg = 2

    # 1. Getting element type
    element_type = defineElemType(elem_type_id)
    if (element_type === nothing)
        @error("Element type passed to fem2D() is unknown")
        return
    end

    # 2. Getting mesh type
    mesh_type = typeMeshFromElement(elem_type_id)

    # Model type
    model_type = default_modeltype

    parameters = ProcessPars(testMaterialProperties(), testBC(), testLoad(), 
        generateTestMesh2D(2), ModelPars(model_type))

    # 3. Reading mesh
    parameters.mesh = read_mesh_from_med(mesh_path, mesh_type)
    
    # 4. Reading problem data
    read_params_JSON!(data_path, parameters)

    if (!isnonloc(parameters))
        @error("Input non-local parameter are not found in elasmech_2d_nonloc()")
        return
    end

    beta_loc = parameters.model.nlpars.locimpact_coeff
    beta_nonloc = parameters.model.nlpars.nonlocimpact_coeff
    impactdist = parameters.model.nlpars.impactdist

    # 5. Integration order
    int_order = 2

    # 6. Problem constants
    nu = parameters.materialProperties[poisC]
    E = parameters.materialProperties[youngMod]
    C = elasticityMatrix(E, nu, parameters.model.type)

    # 7. Local part of stiffness matrix (left part of final equation)
    ensemble_matrix = zeros(2 * size(parameters.mesh.nodes)[1], 2 * 
        size(parameters.mesh.nodes)[1])
    for element_num in eachindex(parameters.mesh.elements)
        k = stiffnessMatrix(C, parameters, element_num, int_order, element_type)
        k .*= beta_loc;
        assembly_left_part!(parameters, ensemble_matrix, k, element_num, freedom_deg)
    end

    # 8. Non-local part of stiffness matrix (left part of final equation)
    for elem_source in eachindex(parameters.mesh.elements)
        # Get list of neighbours
        nodes = parameters.mesh.elements[elem_source]
        neighbours = []
        nodes_cnt = length(nodes)
        x_source = [parameters.mesh.nodes[parameters.mesh.elements[elem_source][i]][1] 
            for i in 1:nodes_cnt]
        y_source = [parameters.mesh.nodes[parameters.mesh.elements[elem_source][i]][2] 
            for i in 1:nodes_cnt]
        startpt_loc = (0, 0)
        startpt_glob = conv_loc_to_glob(startpt_loc[1], startpt_loc[2], x_source, y_source,
            element_type)
        startpt_glob_3d = (startpt_glob[1], startpt_glob[2], 0)
        get_elem_neighbours!(neighbours, elem_source, impactdist, startpt_glob_3d, 
            parameters, element_type)

        # Contribute neighbours impact
        for elem_impact in neighbours
            nonloc_matr = stiffnessmatr_2d_nonloc(C, parameters, elem_source, elem_impact, 
                impactdist, int_order, element_type)
            nonloc_matr .*= beta_nonloc
            # contribute_leftpart_nonloc!(parameters, ensemble_matrix, nonloc_matr, 
                # elem_source, elem_impact, freedom_deg)

            source_connmatr = get_connmatr(parameters, elem_source, freedom_deg)
            impact_connmatr = get_connmatr(parameters, elem_impact, freedom_deg)
            contribute_leftpart_by_connmatr_nonloc!(ensemble_matrix, nonloc_matr, 
                source_connmatr, impact_connmatr)
        end
    end

    # 9. Load vector (right part og final equation)
    load_vector = assembly_loads!(parameters, int_order, element_type, freedom_deg)

    # 10. Applying constraints to equation
    applyConstraints(parameters, load_vector, ensemble_matrix)

    # 11. Creating folder for the equation entities
    mkpath("equation")

    # 12. Writing left part to file
    open("equation/K", "w") do file
        writedlm(file, ensemble_matrix)
    end

    # 13. Writing right part to file
    open("equation/F", "w") do file
        writedlm(file, load_vector)
    end

    # 14. Solving equation
    @info("Solving...")
    @time begin
    result = solve(ensemble_matrix, load_vector)
    end  # @time

    # 15. Exporting result to VTK
    exportToVTK(result, nothing, nothing, nothing, parameters, mesh_type)
end

"""
	elasmech_3d_nonloc(mesh_path::String, data_path::String, elem_type_id::FETypes)

Start calculation of 3D elastic non-local mechanical problem.

# Arguments
- `mesh_path::String`: path to mesh;
- `data_path::String`: path to initial data;
- `elem_type_id::FETypes`: ID of finite element type.
"""
function elasmech_3d_nonloc(mesh_path::String, data_path::String, elem_type_id::FETypes
)
    freedom_deg = 3

    # 1. Getting element type
    element_type = defineElemType(elem_type_id)
    if (element_type === nothing)
        @error("Element type passed to elasmech_3d_nonloc() is unknown")
        return
    end

    # 2. Getting mesh type
    mesh_type = typeMeshFromElement(elem_type_id)

    # Model type
    model_type = default_modeltype

    parameters = ProcessPars(testMaterialProperties(), testBC3D(), testLoad3D(), 
        generateTestMesh3D(), ModelPars(model_type))

    # 3. Reading mesh
    parameters.mesh = read_mesh_from_med(mesh_path, mesh_type)
    
    # 4. Reading problem data
    read_params_JSON!(data_path, parameters)

    if (!isnonloc(parameters))
        @error("Input non-local parameter are not found in elasmech_3d_nonloc()")
        return
    end

    beta_loc = parameters.model.nlpars.locimpact_coeff
    beta_nonloc = parameters.model.nlpars.nonlocimpact_coeff
    impactdist = parameters.model.nlpars.impactdist

    # 5. Integration order
    int_order = 2

    # 6. Problem constants
    nu = parameters.materialProperties[poisC]
    E = parameters.materialProperties[youngMod]
    C = elasticityMatrix(E, nu, parameters.model.type)

    # 7. Local part of stiffness matrix (left part of final equation)
    ensemble_matrix = zeros(3 * size(parameters.mesh.nodes)[1], 3 * 
        size(parameters.mesh.nodes)[1])
    @threads for element_num in eachindex(parameters.mesh.elements)
        k = stiffnessMatrix3D(C, parameters, element_num, int_order, element_type)
        k .*= beta_loc;
        assembly_left_part!(parameters, ensemble_matrix, k, element_num, freedom_deg)
    end

    @info("Assembling non-local part")
    @time begin
    # 8. Non-local part of stiffness matrix (left part of final equation)
    global_neighbours = Dict()
    @info "Searching for neighbours"
    elems_processed = 0
    @threads for elem_source in eachindex(parameters.mesh.elements)
        # Get list of neighbours
        nodes = parameters.mesh.elements[elem_source]
        neighbours = []
        nodes_cnt = length(nodes)
        x_source = [parameters.mesh.nodes[parameters.mesh.elements[elem_source][i]][1] 
            for i in 1:nodes_cnt]
        y_source = [parameters.mesh.nodes[parameters.mesh.elements[elem_source][i]][2] 
            for i in 1:nodes_cnt]
        z_source = [parameters.mesh.nodes[parameters.mesh.elements[elem_source][i]][3] 
            for i in 1:nodes_cnt]
        startpt_loc = (0, 0, 0)
        startpt_glob = conv_loc_to_glob(startpt_loc[1], startpt_loc[2], startpt_loc[3], 
            x_source, y_source, z_source, element_type)
        get_elem_neighbours!(neighbours, elem_source, impactdist, startpt_glob, parameters,
            element_type)

        if isempty(neighbours)
            append!(neighbours, elem_source)
        end

        global_neighbours[elem_source] = neighbours
        elems_processed += 1
        if elems_processed % 10 == 0
            println(elems_processed, " of ", length(parameters.mesh.elements))
        end
    end

    @info "Synchronizing neighbours"
    # Synchronize neighbours
    for (elem, neighbours) in global_neighbours
        for neighbour in neighbours
            counter_neighbours = global_neighbours[neighbour]

            if !(elem in counter_neighbours)
                append!(counter_neighbours, elem)
            end
        end
    end

    @info "Contributing non-local impacts"
    elems_processed = 0
    for (elem_source, neighbours) in global_neighbours
        # Contribute neighbours impact
        @threads for elem_impact in neighbours
            nonloc_matr = stiffnessmatr_3d_nonloc(C, parameters, elem_source, elem_impact, 
                impactdist, int_order, element_type)
            nonloc_matr .*= beta_nonloc
            contribute_leftpart_nonloc!(parameters, ensemble_matrix, nonloc_matr, 
                elem_source, elem_impact, freedom_deg)
        end
        elems_processed += 1
        if elems_processed % 10 == 0
            println(elems_processed, " of ", length(global_neighbours))
        end
    end
    end

    # 9. Load vector (right part og final equation)
    load_vector = assembly_loads!(parameters, int_order, element_type, freedom_deg)

    # 10. Applying constraints to equation
    applyConstraints3D(parameters, load_vector, ensemble_matrix)

    # 11. Creating folder for the equation entities
    # mkpath("equation")

    # 12. Writing left part to file
    # open("equation/K", "w") do file
    #     writedlm(file, ensemble_matrix)
    # end

    # 13. Writing right part to file
    # open("equation/F", "w") do file
    #     writedlm(file, load_vector)
    # end

    # 14. Solving equation
    @info("Solving...")
    @time begin
    result = solve(ensemble_matrix, load_vector)
    end  # @time

    # 15. Verifying result
    # if TestFEM.verify_example(mesh_path, data_path, result)
    #     @info "Result is correct"
    # else
    #     @info "Result is INcorrect"
    # end

    # 16. Writing result to file
    # open("equation/result", "w") do file
    #     writedlm(file, result)
    # end

    # 17. Calclualting deformations
    deformations = calculate_deformations_3d(result, parameters, element_type)

    # 18. Calcualting stresses
    # stresses = calulate_stresses_3d_nonloc(deformations, result, C, beta_loc, beta_nonloc,
        # global_neighbours, impactdist, parameters, int_order, element_type)

    # 19. Calculating Von-Mises stresses
    # von_mises = calculateVonMises(stresses)

    # 20. Exporting result to VTK
    exportToVTK(result, deformations, nothing, nothing, parameters, mesh_type)
end

end  # Core