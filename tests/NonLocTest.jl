using CoreFEM
using MeshFEM

"""
    get_elem_neighbours_test()

Test `CoreFEM.get_elem_neighbours!()` function.
"""
##
function get_elem_neighbours_test()
    parameters = ProcessPars(testMaterialProperties(), testBC(), testLoad(), 
        MeshFEM.generateTestMesh2D(5))
    @info("CoreFEM.get_elem_neighbours!() test...")

    neighbours = [12]
    startpt = (50, 50, 0)

    # Test 1
    CoreFEM.get_elem_neighbours!(neighbours, 12, 10, startpt, parameters)
    print("Test1............")
    if length(neighbours) == 1 && neighbours[1] == 12
        println("Passed")
    else
        println("Failed")
    end

    # Test2
    CoreFEM.get_elem_neighbours!(neighbours, 12, 20, startpt, parameters)
    print("Test2............")
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

"""
    nonloc_gaussimpact()

Test `CoreFEM.nonloc_gaussimpact()` function.
"""
##
function nonloc_gaussimpact()
    @info("CoreFEM.nonloc_gaussimpact() test...")

    # Test 1
    distance = 1
    impactdistance = 1
    normfactor = 1 / (2 * pi * impactdistance^2)
    impact = CoreFEM.nonloc_gaussimpact(normfactor, impactdistance, distance)
    print("Test1............")
    if impact ≈ (1 / (2 * ℯ * pi))
        println("Passed")
    else
        println("Failed")
    end

    # Test 2
    distance = 1
    impactdistance = 2
    normfactor = 1 / (2 * pi * impactdistance^2)
    impact = CoreFEM.nonloc_gaussimpact(normfactor, impactdistance, distance)
    print("Test2............")
    if impact ≈ (1 / (8 * exp(1 / 4) * pi))
        println("Passed")
    else
        println("Failed")
    end

    # Test 3
    distance = 2
    impactdistance = 6
    normfactor = 1 / (2 * pi * impactdistance^2)
    impact = CoreFEM.nonloc_gaussimpact(normfactor, impactdistance, distance)
    print("Test3............")
    if impact ≈ (1 / (72 * exp(1 / 9) * pi))
        println("Passed")
    else
        println("Failed")
    end
end

nonloc_gaussimpact()
##
