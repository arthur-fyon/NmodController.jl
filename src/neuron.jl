# Prebuild gating functions
boltz(V, A, B) = 1 / (1 + exp((V+A)/B))
tauX(V, A, B, D, E) = A - B / (1+exp((V+D)/E))

# Data structure that defines a neuronal current
struct IonCurrent
    name::String
    numberOfGatings::Int64
    exponents::Vector{Int64}
    steadyStateGatings::Vector{Function}
    timeConstants::Vector{Any}
    reversalPotential::Float64
    calciumDependency::Bool
    MgDependency::Bool
end

# Data structure that defines calcium calcium dynamics
struct CalciumDynamic
    numberOfCurrents::Int64
    currentNames::Vector{String}
    coefficients::Vector{Float64}
    nernstEquilibriumValue::Float64
    timeConstant::Float64
end

# Data structure that defines the neuronal model
struct NeuronCB
    C::Float64
    numberOfIonCurrents::Int64
    ionCurrents::Vector{IonCurrent}
    maximumConductances::Vector{Float64}
    leakageConductance::Float64
    reversaleLeakagePotential::Float64
    calciumDynamics::Union{CalciumDynamic, Bool}
    globalCalciumDependency::Bool
    globalMgDependency::Bool
    Vth::Float64
end

# Function that initializes a certain type of current
function initializeCurrent(name::String, reversalPotential::Union{Int64, Float64}; numberOfGatings::Int64=1, exponents::Union{Int64, Vector{Int64}}=1,
    activationSteadyStateGating::Union{Function, Float64}=0., activationTimeConstants::Union{Function, Float64, Int64}=0.,
    inactivationSteadyStateGating::Union{Function, Float64}=0., inactivationTimeConstants::Union{Function, Float64, Int64}=0.,
    calciumDependency::Bool=false, MgDependency::Bool=false)

    # If no activation steady state function is provided, throw error
    if !(typeof(activationSteadyStateGating) <: Function)
        error("Please enter a steady state activation gating value of type Function!")
    end

    # If we only have an activation gating variable
    if numberOfGatings == 1
        # If more than one exponent is provided, throw error
        if !isa(exponents, Int64)
            error("Given that you only have one gating variable in the ionic current, exponents must be of type Int64!")
        end

        # Create the IonCurrent structure
        return IonCurrent(name, numberOfGatings, [exponents, ], [activationSteadyStateGating, ],
            [activationTimeConstants, ], Float64(reversalPotential), calciumDependency, MgDependency)

    # If we have both an activation and inactivation gating variables
    elseif numberOfGatings == 2
        # If anything else than 2 exponents are provided, throw error
        if !isa(exponents, Vector{Int64}) && length(exponents) ≠ 2
            error("Given that you have two gating variables in the ionic current, exponents must be of type Vector{Int64}!")
        end
        # If no inactivation steady state function is provided, throw error
        if !(typeof(inactivationSteadyStateGating) <: Function)
            error("Please provide a steady state inactivation gating value of type Function!")
        end

        # Create the IonCurrent structure
        return IonCurrent(name, numberOfGatings, exponents, [activationSteadyStateGating, inactivationSteadyStateGating],
            [activationTimeConstants, inactivationTimeConstants], Float64(reversalPotential), calciumDependency, MgDependency)

    # If less than 1 or more than 2 gating variables are provided, throw error
    else
        error("Number of gating variables must be equal to 1 or 2!")
    end
end

# Function that initializes calcium dynamics
function initializeCalciumDynamics(currentNames::Union{String, Vector{String}}, coefficients::Union{Int64, Float64, Vector{Int64}, Vector{Float64}},
    nernstEquilibriumValue::Union{Int64, Float64}, timeConstant::Union{Int64, Float64})

    # If there is only one current in the calcium dynamics
    if (isa(coefficients, Int64) || isa(coefficients, Float64)) && isa(currentNames, String)
        # Create the CalciumDynamic structure
        return CalciumDynamic(1, [currentNames, ], [Float64(coefficients), ], Float64(nernstEquilibriumValue), Float64(timeConstant))

    # If there is more than one current in the calcium dynamics
    elseif (isa(coefficients, Vector{Int64}) || isa(coefficients, Vector{Float64})) && isa(currentNames, Vector{String})
        # If length of currentNames and coefficients are the same
        if length(coefficients) == length(currentNames)
            # Create the CalciumDynamic structure
            return CalciumDynamic(length(coefficients), currentNames, Vector{Float64}(coefficients), Float64(nernstEquilibriumValue), Float64(timeConstant))
        # If length of currentNames and coefficients are not the same, throw error
        else
            error("Length of currentNames and coefficients must be the same!")
        end

    # Else, throw an error
    else
        error("currentNames and coefficients must be both of type Vector or of type String and Float64/Int64 respectively!")
    end
end

# Function that initializes neuronal model
function initializeNeuronModel(ionCurrents::Union{IonCurrent, Vector{IonCurrent}}; C::Float64=1., leakageConductance::Union{Float64, Int64}=1., 
    reversaleLeakagePotential::Union{Float64, Int64}=-50., calciumDynamics::Union{CalciumDynamic, Bool}=false, 
    maximumConductances::Union{Vector{Int64}, Int64, Vector{Float64}, Float64}=0.)

    # If there is only one current in the model, throw error
    if isa(ionCurrents, IonCurrent)
        error("ionCurrents must be of type Vector{IonCurrent} of minimum length 2, i.e. at least 2 ion currents!")
    elseif length(ionCurrents) <= 1
        error("ionCurrents must be of type Vector{IonCurrent} of minimum length 2, i.e. at least 2 ion currents!")
    end

    # Extracting the number of ion currents in the model
    numberOfIonCurrents = length(ionCurrents)

    # Checking whether there is somewhere calcium or mg dynamics in the model
    globalCalciumDependency = false
    globalMgDependency = false
    for ionCurr in ionCurrents
        if ionCurr.calciumDependency
            globalCalciumDependency = true
        elseif ionCurr.MgDependency
            globalMgDependency = false
        end
    end

    # If there is calcium dependency in the model but no calcium dynamics are provided, throw error
    if globalCalciumDependency && isa(calciumDynamics, Bool)
        error("Please provide a calciumDynamics of type CalciumDynamic!")
    end

    # If the user did not enter any maximumConductances, fill a vector with NaN
    if isa(maximumConductances, Float64) && maximumConductances == 0.
        maximumConductances = zeros(numberOfIonCurrents) .+ NaN
    # If the length of the maximum ion channels is the same as the one of the ion currents structure, throw error
    elseif (isa(maximumConductances, Vector{Float64}) || isa(maximumConductances, Vector{Int64})) && length(maximumConductances) ≠ numberOfIonCurrents
        error("Length of ionCurrents and maximumConductances must be the same!")
    # If the maximumConductances is not the of type Vector{Float64}, throw error
    elseif (isa(maximumConductances, Float64) && maximumConductances ≠ 0.) || isa(maximumConductances, Int64)
        error("maximumConductances must be of type Vector{Float64}!")
    end

    # Create the NeuronCB structure
    return NeuronCB(C, numberOfIonCurrents, ionCurrents, Vector{Float64}(maximumConductances), Float64(leakageConductance), 
        Float64(reversaleLeakagePotential), calciumDynamics, globalCalciumDependency, globalMgDependency, NaN)
end
