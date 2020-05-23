using CoreFEM

# Method assumes that result is for 2D object
"""
    exportToCSV(result::Array, pars::processPars)

Export given result to CSV file.

# Arguments
- `result::Array`: results vector (assuming that given result is for some 2D object);
- `pars::processPars`: parameters of current model.
"""
function exportToCSV(result::Array, pars::processPars)
    displacementHeader = "X, Y, Z, Displacement\n"
    open("output/displacementX.csv", "w") do file
        write(file, displacementHeader)
        for i in 1:size(pars.mesh.nodes)[1]
            write(file, string(pars.mesh.nodes[i][1]) * ", " * string(pars.mesh.nodes[i][2]) * ", " * "0, " * string(result[2 * i - 1]) * "\n")
        end
    end
    open("output/displacementY.csv", "w") do file
        write(file, displacementHeader)
        for i in 1:size(pars.mesh.nodes)[1]
            write(file, string(pars.mesh.nodes[i][1]) * ", " * string(pars.mesh.nodes[i][2]) * ", " * "0, " * string(result[2 * i]) * "\n")
        end
    end
end  # exportToCSV

function exportToDAT2D(result::Array, pars::processPars)
    open("output/displacementX.dat", "w") do file
        for i in 1:size(pars.mesh.nodes)[1]
            write(file, string(pars.mesh.nodes[i][1]) * " " * string(pars.mesh.nodes[i][2]) * " " * string(result[2 * i - 1]) * "\n")
        end
    end
    open("output/displacementY.dat", "w") do file
        for i in 1:size(pars.mesh.nodes)[1]
            write(file, string(pars.mesh.nodes[i][1]) * " " * string(pars.mesh.nodes[i][2]) * " " * string(result[2 * i]) * "\n")
        end
    end
end # exportToDAT2D