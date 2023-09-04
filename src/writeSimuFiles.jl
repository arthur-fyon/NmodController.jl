# Function that writes conductance-based model ODE equations in a jl files
function writeUncontrolledODEs(neuron::NeuronCB)
    # First removing the file if it already exists
    if isfile("CBModelODEs.jl")
        rm("CBModelODEs.jl")
    end

    # Opening file in write mode
    f = open("CBModelODEs.jl", "w")

    # Writing the first commentary line
    line = "#=\n\t"
    write(f, line)
    line = "This file contains differential equations describing the CB model of interest\n"
    write(f, line)
    line = "=#"
    write(f, line)

    write(f, "\n")

    # Writing the function prototype
    line = "# Function that outputs values of variables derivatives\n"
    write(f, line)
    line = "function CB_ODE(dx, x, p, t)\n\t"
    write(f, line)

    # Writing parameters naming
    line = "# Parameters\n\t"
    write(f, line)
    line = "Iapp = p[1](t) # Time dependent applied current\n\t"
    write(f, line)
    
    # Looping over all model currents
    for (i, ionCurrent) in enumerate(neuron.ionCurrents)
        line = "g"*ionCurrent.name*" = p["*string(i+1)*"] # Maximum conductance of current "*ionCurrent.name*"\n\t"
        write(f, line)
    end

    # Leakage conductance and capacitance
    line = "gleak = p["*string(length(neuron.ionCurrents)+2)*"] # Leakage conductance\n\t"
    write(f, line)
    line = "C = p["*string(length(neuron.ionCurrents)+3)*"] # Membrane capacitance\n\t"
    write(f, line)

    # If there is Mg dependency, add Mg to parameter list
    if neuron.globalMgDependency
        line = "Mg = p["*string(length(neuron.ionCurrents)+4)*"] # Mg concentration\n\t"
        write(f, line)
    end

    write(f, "\n\t")

    # Writing variables naming
    line = "# Variables\n\t"
    write(f, line)
    line = "V = x[1](t)\n\t"
    write(f, line)

    # Looping over all model currents
    i = 2
    for ionCurrent in neuron.ionCurrents
        # Writing the activation variable
        line = "m"*ionCurrent.name*" = x["*string(i)*"] # Activation gating variable of current "*ionCurrent.name*"\n\t"
        write(f, line)
        i = i + 1

        # Add an inactivation variable if needed
        if ionCurrent.numberOfGatings > 1 
            line = "h"*ionCurrent.name*" = x["*string(i)*"] # Inactivation gating variable of current "*ionCurrent.name*"\n\t"
            write(f, line)
            i = i + 1
        end
    end

    # If there is Ca dependency, add Ca to variable list
    if neuron.globalCalciumDependency
        line = "Ca = x["*string(i)*"] # Ca concentration\n\t"
        write(f, line)
    end

    write(f, "\n\t")

    # Writing ODEs
    line = "# ODEs\n\t"
    write(f, line)

    # Writing voltage ODE
    line = "dx[1] = (1/C) * (- "
    write(f, line)

    # Looping over all model currents
    for ionCurrent in neuron.ionCurrents
        if ionCurrent.numberOfGatings == 1 
            line = "g"*ionCurrent.name*"*m"*ionCurrent.name*"^"*string(ionCurrent.exponents[1])*"*(V - "*string(ionCurrent.reversalPotential)*") - \n\t\t\t\t\t   "
            write(f, line)
        else
            line = "g"*ionCurrent.name*"*m"*ionCurrent.name*"^"*string(ionCurrent.exponents[1])*"*h"*ionCurrent.name*"^"*string(ionCurrent.exponents[2])*
                "*(V - "*string(ionCurrent.reversalPotential)*") - \n\t\t\t\t\t   "
            write(f, line)
        end
    end

    # Leakage conductance and applied current
    line = "gleak*(V - "*string(neuron.reversaleLeakagePotential)*") + Iapp)\n\t"
    write(f, line)

    # Writing gating variable ODEs by looping over all model current
    i = 2
    for ionCurrent in neuron.ionCurrents
        if ionCurrent.calciumDependency
            # Writing the activation variable
            line = "dx["*string(i)*"] = (1/tau_m"*ionCurrent.name*"(V)) * (m"*ionCurrent.name*"_inf(V, Ca) - m"*ionCurrent.name*")\n\t"
            write(f, line)
            i = i + 1

            # Add an inactivation variable if needed
            if ionCurrent.numberOfGatings > 1 
                line = "dx["*string(i)*"] = (1/tau_h"*ionCurrent.name*"(V)) * (h"*ionCurrent.name*"_inf(V, Ca) - h"*ionCurrent.name*")\n\t"
                write(f, line)
                i = i + 1
            end
        elseif ionCurrent.MgDependency
            # Writing the activation variable
            line = "dx["*string(i)*"] = (1/tau_m"*ionCurrent.name*"(V)) * (m"*ionCurrent.name*"_inf(V, Mg) - m"*ionCurrent.name*")\n\t"
            write(f, line)
            i = i + 1

            # Add an inactivation variable if needed
            if ionCurrent.numberOfGatings > 1 
                line = "dx["*string(i)*"] = (1/tau_h"*ionCurrent.name*"(V)) * (h"*ionCurrent.name*"_inf(V, Mg) - h"*ionCurrent.name*")\n\t"
                write(f, line)
                i = i + 1
            end
        else
            # Writing the activation variable
            line = "dx["*string(i)*"] = (1/tau_m"*ionCurrent.name*"(V)) * (m"*ionCurrent.name*"_inf(V) - m"*ionCurrent.name*")\n\t"
            write(f, line)
            i = i + 1

            # Add an inactivation variable if needed
            if ionCurrent.numberOfGatings > 1 
                line = "dx["*string(i)*"] = (1/tau_h"*ionCurrent.name*"(V)) * (h"*ionCurrent.name*"_inf(V) - h"*ionCurrent.name*")\n\t"
                write(f, line)
                i = i + 1
            end
        end
    end

    # Writing calcium dynamics if needed
    if neuron.globalCalciumDependency
        # Extracting calcium dynamics
        calciumDynamics = neuron.calciumDynamics

        # Extracting the current of interest
        if length(calciumDynamics.coefficients) == 1
            currentIndex = findfirst(x -> occursin(calciumDynamics.currentNames[1], x.name), neuron.ionCurrents)
            if isa(currentIndex, Nothing)
                error("Current of the calcium dynamics not found in the model!")
            end
            ionCurrentCa = [neuron.ionCurrents[currentIndex], ]
        elseif length(calciumDynamics.coefficients) == 2
            currentIndices = zeros(Int64, 2)
            for (i, currentName) in enumerate(calciumDynamics.currentNames)
                currentIndex = findfirst(x -> occursin(currentName, x.name), neuron.ionCurrents)
                if isa(currentIndex, Nothing)
                    error("Current of the calcium dynamics not found in the model!")
                end
                currentIndices[i] = currentIndex
            end
            if isa(currentIndices, Nothing)
                error("Current of the calcium dynamics not found in the model!")
            end

            ionCurrentCa = [neuron.ionCurrents[currentIndices[1]], neuron.ionCurrents[currentIndices[2]]]
        end

        # Writing by looping over all Ca currents
        line = "dx["*string(i)*"] = ("
        write(f, line)
        
        # Looping over all Ca currents of the model
        ionCurrentCa = neuron.ionCurrents[currentIndices]
        for (i, ionCurrent) in enumerate(ionCurrentCa)
            line = string(calciumDynamics.coefficients[i])
            write(f, line)

            if ionCurrent.numberOfGatings == 1 
                line = "*g"*ionCurrent.name*"*m"*ionCurrent.name*"^"*string(ionCurrent.exponents[1])*"*(V - "*string(ionCurrent.reversalPotential)*") + "
                write(f, line)
            else
                line = "*g"*ionCurrent.name*"*m"*ionCurrent.name*"^"*string(ionCurrent.exponents[1])*"*h"*ionCurrent.name*"^"*string(ionCurrent.exponents[2])*
                    "*(V - "*string(ionCurrent.reversalPotential)*") + "
                write(f, line)
            end
        end

        # Adding the nernst potential Ca value
        line = "-Ca + "*string(calciumDynamics.nernstEquilibriumValue)*") / "*string(calciumDynamics.timeConstant)*"\n"
        write(f, line)
    end

    # Finishing the ODE function
    line = "end\n"
    write(f, line)

    # Closing the file
    close(f)
end