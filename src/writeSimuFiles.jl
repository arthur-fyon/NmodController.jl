"""
    writeUncontrolledODEs(neuron::NeuronCB; 
        filename::String="CBModelODEs.jl")

Write a julia file that contains all ODEs for a specific conductance based model. The equations can be integrated using the DifferentialEquations.jl package.
Returns nothing.

# Arguments
- `neuron`: data structure containing the conductance based model of interest.
- `filename`: name of the julia file that will contain the equations. Optional.

# Example
```jldoctest
julia> writeUncontrolledODEs(STG)
```
"""
function writeUncontrolledODEs(neuron::NeuronCB; filename::String="CBModelODEs.jl")
    # If filename last 3 carachters are not ".jl", throw error
    if filename[end-2:end] ≠ ".jl"
        error("Filename must be a Julia file, i.e. being terminated by .jl!")
    end

    # First removing the file if it already exists
    if isfile(filename)
        rm(filename)
    end

    # Opening file in write mode
    f = open(filename, "w")

    # Writing the first commentary line
    line = "#=\n\t"
    write(f, line)
    line = "This file contains differential equations describing the CB model of interest\n"
    write(f, line)
    line = "=#\n"
    write(f, line)

    write(f, "\n")

    # Writing the function prototype
    line = "# Function that outputs values of variables derivatives\n"
    write(f, line)
    line = "function "*filename[1:end-3]*"(dx, x, p, t)\n\t"
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
    line = "V = x[1] # Membrane voltage\n\t"
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
            for (j, currentName) in enumerate(calciumDynamics.currentNames)
                currentIndex = findfirst(x -> occursin(currentName, x.name), neuron.ionCurrents)
                if isa(currentIndex, Nothing)
                    error("Current of the calcium dynamics not found in the model!")
                end
                currentIndices[j] = currentIndex
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
        for (j, ionCurrent) in enumerate(ionCurrentCa)
            line = string(calciumDynamics.coefficients[j])
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

"""
    writeControlledODEs(neuron::NeuronCB, 
        controlledConductances::Vector{String}, 
        controlledDICs::Vector{String}; 
        filename::String="ControlledCBModelODEs.jl")

Write a julia file that contains all ODEs for a specific conductance based model that is coupled with a neuromodulation controller. The equations can be integrated using the DifferentialEquations.jl package.
Refer to Fyon et al., 2023 "Reliable neuromodulation from adaptive control of ion channel expression". Returns nothing.

# Arguments
- `neuron`: data structure containing the conductance based model of interest.
- `controlledConductances`: vector of strings that contains the names of the conductances the controller outputs.
- `controlledDICs`: vector of strings that contains the names of the DICs the controller gets in input.
- `filename`: name of the julia file that will contain the equations. Optional.

# Example
```jldoctest
julia> writeControlledODEs(STG, ["CaS", "A"], ["s", "u"])
```
"""
function writeControlledODEs(neuron::NeuronCB, controlledConductances::Vector{String}, controlledDICs::Vector{String}; filename::String="ControlledCBModelODEs.jl")
    # If there are more DICs controlled than modulated conductances, throw error
    if length(controlledDICs) > length(controlledConductances)
        error("There should be less controlled neuronal feedback gains than modulated conductances!")
    end

    # If filename last 3 carachters are not ".jl", throw error
    if filename[end-2:end] ≠ ".jl"
        error("Filename must be a Julia file, i.e. being terminated by .jl!")
    end

    # First removing the file if it already exists
    if isfile(filename)
        rm(filename)
    end

    # First retrieving the indices of the controlled conductances
    controlledIndices = []
    for (i, ionCurrent) in enumerate(neuron.ionCurrents)
        for controlledConductance in controlledConductances
            if controlledConductance == ionCurrent.name
                push!(controlledIndices, i)
            end
        end
    end

    # If no controlled conductances in the model, throw error
    if isa(controlledIndices, Nothing)
        error("No controlled conductances found in the model!")
    end

    # Opening file in write mode
    f = open(filename, "w")

    # Writing the first commentary line
    line = "#=\n\t"
    write(f, line)
    line = "This file contains differential equations describing the controlled CB model of interest\n"
    write(f, line)
    line = "=#\n"
    write(f, line)

    write(f, "\n")

    # Writing the function prototype
    line = "# Function that outputs values of variables derivatives\n"
    write(f, line)
    line = "function "*filename[1:end-3]*"(dx, x, p, t)\n\t"
    write(f, line)

    # Writing parameters naming
    line = "# Parameters\n\t"
    write(f, line)
    line = "Iapp = p[1](t) # Time dependent applied current\n\t"
    write(f, line)
    
    # Looping over all model currents
    i = 2
    for (j, ionCurrent) in enumerate(neuron.ionCurrents)
        if !(j in controlledIndices)
            line = "g"*ionCurrent.name*" = p["*string(i)*"] # Maximum conductance of current "*ionCurrent.name*"\n\t"
            i = i + 1
            write(f, line)
        end
    end

    # Leakage conductance and capacitance
    line = "gleak = p["*string(i)*"] # Leakage conductance\n\t"
    i = i + 1
    write(f, line)
    line = "C = p["*string(i)*"] # Membrane capacitance\n\t"
    i = i + 1
    write(f, line)

    # If there is Mg dependency, add Mg to parameter list
    if neuron.globalMgDependency
        line = "Mg = p["*string(i)*"] # Mg concentration\n\t"
        i = i + 1
        write(f, line)
    end

    # Controller parameters
    line = "α = p["*string(i)*"] # Rate of transfer between intracellular and membrane\n\t"
    i = i + 1
    write(f, line)
    line = "β = p["*string(i)*"] # Rate of degradation of intracellular proteins\n\t"
    i = i + 1
    write(f, line)
    line = "Kp = p["*string(i)*"] # Proportional gain\n\t"
    i = i + 1
    write(f, line)
    line = "Ki = p["*string(i)*"] # Integral gain\n\t"
    i = i + 1
    write(f, line)
    line = "SVth = p["*string(i)*"] # Sensitivity matrix at threshold voltage\n\t"
    i = i + 1
    write(f, line)

    # DICs reference
    if "f" in controlledDICs
        line = "gfth = p["*string(i)*"](t) # Reference gf(Vth)\n\t"
        i = i + 1
        write(f, line)
    end
    if "s" in controlledDICs
        line = "gsth = p["*string(i)*"](t) # Reference gs(Vth)\n\t"
        i = i + 1
        write(f, line)
    end
    if "u" in controlledDICs
        line = "guth = p["*string(i)*"](t) # Reference gu(Vth)\n\t"
        i = i + 1
        write(f, line)
    end

    write(f, "\n\t")

    # Writing variables naming
    line = "# Variables\n\t"
    write(f, line)
    line = "V = x[1] # Membrane voltage\n\t"
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
        i = i + 1
    end

    # Variables of the controller
    for controlledCurrent in neuron.ionCurrents[controlledIndices]
        line = "g"*controlledCurrent.name*"i = x["*string(i)*"] # Intracellular maximum conductance of current "*controlledCurrent.name*"\n\t"
        write(f, line)
        i = i + 1

        line = "g"*controlledCurrent.name*" = x["*string(i)*"] # Maximum conductance of current "*controlledCurrent.name*"\n\t"
        write(f, line)
        i = i + 1

        line = "z"*controlledCurrent.name*" = x["*string(i)*"] # Integral variable of current "*controlledCurrent.name*"\n\t"
        write(f, line)
        i = i + 1
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
            for (j, currentName) in enumerate(calciumDynamics.currentNames)
                currentIndex = findfirst(x -> occursin(currentName, x.name), neuron.ionCurrents)
                if isa(currentIndex, Nothing)
                    error("Current of the calcium dynamics not found in the model!")
                end
                currentIndices[j] = currentIndex
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
        for (j, ionCurrent) in enumerate(ionCurrentCa)
            line = string(calciumDynamics.coefficients[j])
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
        line = "-Ca + "*string(calciumDynamics.nernstEquilibriumValue)*") / "*string(calciumDynamics.timeConstant)*"\n\t"
        i = i + 1
        write(f, line)
    end

    # Retrieving which lines of the sensitivty matrix we need
    timescales = []
    if "f" in controlledDICs
        push!(timescales, 1)
    end

    if "s" in controlledDICs
        push!(timescales, 2)
    end

    if "u" in controlledDICs
        push!(timescales, 3)
    end

    write(f, "\n\t")
    line = "# Retrieving which line of the sensitivity matrix matter\n\t"
    write(f, line)

    if length(timescales) == 1
        line = "timescales = ["*string(timescales[1])*"]\n\t"
        write(f, line)
    else
        line = "timescales = ["
        write(f, line)
        for timescale in timescales[1:end-1]
            line = string(timescale)*", "
            write(f, line)
        end
        line = string(timescales[end])*"]\n\t"
        write(f, line)
    end

    # Retrieving which column belong to unmodulated conductances
    write(f, "\n\t")
    line = "# Retrieving which column of the sensitivity matrix belong to unmodulated conductances\n\t"
    write(f, line)
    unmodulated = Vector(1 : length(neuron.ionCurrents))
    for controlledIndex in controlledIndices
        deleteat!(unmodulated, findall(x -> x == controlledIndex, unmodulated))
    end

    if length(unmodulated) == 1
        line = "unmodulated = ["*string(unmodulated[1])*"]\n\t"
        write(f, line)
    else
        line = "unmodulated = ["
        write(f, line)
        for unmodulated_i in unmodulated[1:end-1]
            line = string(unmodulated_i)*", "
            write(f, line)
        end
        line = string(unmodulated[end])*"]\n\t"
        write(f, line)
    end

    # Retrieving which column belong to modulated conductances
    write(f, "\n\t")
    line = "# Retrieving which column of the sensitivity matrix belong to modulated conductances\n\t"
    write(f, line)
    if length(controlledIndices) == 1
        line = "modulated = ["*string(controlledIndices[1])*"]\n\t"
        write(f, line)
    else
        line = "modulated = ["
        write(f, line)
        for controlledIndex in controlledIndices[1:end-1]
            line = string(controlledIndex)*", "
            write(f, line)
        end
        line = string(controlledIndices[end])*"]\n\t"
        write(f, line)
    end

    # Computing the right hand side of the linear system
    write(f, "\n\t")
    line = "# Computing the right hand side of the linear system\n\t"
    write(f, line)
    if length(timescales) == 1
        if 1 in timescales
            line = "gDICr = [gfth, ]\n\t"
            write(f, line)
        end
        if 2 in timescales
            line = "gDICr = [gsth, ]\n\t"
            write(f, line)
        end
        if 3 in timescales
            line = "gDICr = [guth, ]\n\t"
            write(f, line)
        end
        
    else
        line = "gDICr = ["
        write(f, line)
        for timescale in timescales[1:end-1]
            if timescale == 1
                line = "gfth, "
                write(f, line)
            elseif timescale == 2
                line = "gsth, "
                write(f, line)
            elseif timescale == 2
                line = "guth, "
                write(f, line)
            end
        end
        if timescales[end] == 1
            line = "gfth]\n\t"
            write(f, line)
        elseif timescales[end] == 2
            line = "gsth]\n\t"
            write(f, line)
        elseif timescales[end] == 3
            line = "guth]\n\t"
            write(f, line)
        end
    end

    line = "gDICr = gDICr - SVth[timescales, unmodulated] * collect(p[2:"*string(length(neuron.ionCurrents)-length(controlledConductances)+1)*"])\n\t"
    write(f, line)

    # Computing the left hand side of the linear system
    write(f, "\n\t")
    line = "# Computing the left hand side of the linear system\n\t"
    write(f, line)
    line = "Smod = SVth[timescales, modulated]\n\t"
    write(f, line)

    write(f, "\n\t")
    line = "# Computing the solution of the linear system\n\t"
    write(f, line)
    if length(timescales) == length(controlledIndices)
        line = "g_r = \\(Smod, gDICr)\n\t"
        write(f, line)
    else
        line = "g_r = transpose(Smod) * inv(Smod * transpose(Smod)) * gDICr\n\t"
        write(f, line)
    end

    # Computing signal errors and control inputs
    write(f, "\n\t")
    line = "# Error signals and control inputs\n\t"
    write(f, line)
    for (j, controlledIndex) in enumerate(controlledIndices)
        line = "e"*neuron.ionCurrents[controlledIndex].name*" = g_r["*string(j)*"] - g"*neuron.ionCurrents[controlledIndex].name*"\n\t"
        write(f, line)

        line = "u"*neuron.ionCurrents[controlledIndex].name*" = Kp * e"*neuron.ionCurrents[controlledIndex].name*" + Ki * z"*neuron.ionCurrents[controlledIndex].name*"\n\t"
        write(f, line)
    end

    # ODEs of the controller
    write(f, "\n\t")
    line = "# ODEs of the controller\n\t"
    write(f, line)
    for controlledIndex in controlledIndices
        line = "dx["*string(i)*"] = α * g"*neuron.ionCurrents[controlledIndex].name*" - α * g"*neuron.ionCurrents[controlledIndex].name*"i - β * g"*neuron.ionCurrents[controlledIndex].name*"i + u"*neuron.ionCurrents[controlledIndex].name*"\n\t"
        i = i + 1
        write(f, line)

        line = "dx["*string(i)*"] = α * g"*neuron.ionCurrents[controlledIndex].name*"i - α * g"*neuron.ionCurrents[controlledIndex].name*"\n\t"
        i = i + 1
        write(f, line)

        line = "dx["*string(i)*"] = e"*neuron.ionCurrents[controlledIndex].name*"\n\t"
        i = i + 1
        write(f, line)
    end
    write(f, "\n")

    # Finishing the ODE function
    line = "end\n"
    write(f, line)

    # Closing the file
    close(f)
end
