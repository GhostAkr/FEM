using CoreFEM
using MeshFEM

"""
    get_elem_neighbours_test()

Test `CoreFEM.get_elem_neighbours!()` function.
"""
##
function get_elem_neighbours_test()
    parameters = processPars(testMaterialProperties(), testBC(), testLoad(), 
        MeshFEM.generateTestMesh2D(5))
    @info("CoreFEM.get_elem_neighbours!() test...")

    neighbours = [12]
    startpt = (50, 50, 0)

    # Test 1
    CoreFEM.get_elem_neighbours!(neighbours, 12, 10, startpt, parameters)
    print("Test............")
    if length(neighbours) == 1 && neighbours[1] == 12
        println("Passed")
    else
        println("Failed")
    end

    # Test2
    CoreFEM.get_elem_neighbours!(neighbours, 12, 20, startpt, parameters)
    print("Test............")
    if length(neighbours) == 9 &&
            6 in neighbours &&
            11 in neighbours &&
            16 in neighbours &&
            17 in neighbours &&
            18 in neighbours &&
            13 in neighbours &&
            8 in neighbours &&
            7 in neighbours &&
            12 in neighbours
        println("Passed")
    else
        println("Failed")
    end
end

get_elem_neighbours_test()
##
