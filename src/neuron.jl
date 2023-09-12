"
    struct IonCurrent
        name::String
        numberOfGatings::Int64
        exponents::Vector{Int64}
        steadyStateGatings::Vector{Function}
        timeConstants::Vector{Function}
        reversalPotential::Float64
        calciumDependency::Bool
        MgDependency::Bool
    end

Data structure defining an ionic current in a conductance based model.

To initialize a certain type of current, please use the function initializeCurrent.
"
struct IonCurrent
    "Name of the ionic current"
    name::String
    "Number of gating variables in the ionic current"
    numberOfGatings::Int64
    "Exponents of gating variables"
    exponents::Vector{Int64}
    "Equilibrium functions of gating variables"
    steadyStateGatings::Vector{Function}
    "Time constant functions of gating variables"
    timeConstants::Vector{Function}
    "Reversal Nernst potential of the ion"
    reversalPotential::Float64
    "Flag indicating if the ionic current depends on the calcium"
    calciumDependency::Bool
    "Flag indicating if the ionic current depends on the magnesium"
    MgDependency::Bool
end

"
    struct CalciumDynamic
        numberOfCurrents::Int64
        currentNames::Vector{String}
        coefficients::Vector{Float64}
        nernstEquilibriumValue::Float64
        timeConstant::Float64
    end

Data structure defining an ODE for intracellular calcium dynamics.

To initialize a certain calcium dynamics, please use the function initializeCalciumDynamics.
"
struct CalciumDynamic
    "Number of calcium ionic current involved in the ODE"
    numberOfCurrents::Int64
    "Name(s) of the calcium ionic current(s) involved in the ODE"
    currentNames::Vector{String}
    "Coefficient(s) of the calcium ionic current(s) involved in the ODE"
    coefficients::Vector{Float64}
    "Equilibrium value of the intracellular calcium without any current"
    nernstEquilibriumValue::Float64
    "Time constant of the intracellular calcium"
    timeConstant::Float64
end

"
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
    end

Data structure defining a complete conductance based model.

To initialize a certain type of model, please use the function initializeNeuronModel.
"
struct NeuronCB
    "Capacitance of the membrane"
    C::Float64
    "Number of active ionic current in the model (without the leakage one)"
    numberOfIonCurrents::Int64
    "Vector of the ionic current data structures"
    ionCurrents::Vector{IonCurrent}
    "Vector of the maximum ion channel conductances (without the leakage one)"
    maximumConductances::Vector{Float64}
    "Leakage conductance"
    leakageConductance::Float64
    "Reversal Nernst potential for the leakage current"
    reversaleLeakagePotential::Float64
    "Data structure for describing intracellular calcium dynamics ODE"
    calciumDynamics::Union{CalciumDynamic, Bool}
    "Flag indicating if the model depends on the calcium"
    globalCalciumDependency::Bool
    "Flag indicating if the model depends on the magnesium"
    globalMgDependency::Bool
end

"""
    initializeCurrent(name::String, reversalPotential::Union{Int64, Float64}; numberOfGatings::Int64=1, exponents::Union{Int64, Vector{Int64}}=1,
    activationSteadyStateGating::Union{Function, Float64}=0., activationTimeConstant::Union{Function, Float64, Int64}=NaN,
    inactivationSteadyStateGating::Union{Function, Float64}=0., inactivationTimeConstant::Union{Function, Float64, Int64}=NaN,
    calciumDependency::Bool=false, MgDependency::Bool=false)

Initialize an ionic current data structure. The current must be of type Iion = bar_gion * m^a1 * h^a2 * (V - Eion).

Note that the maximum ion channel conductance bar_gion is not contained in the ionic current data structure but in the neuronal model one.

# Arguments
- `name::String`: name of the ionic current.
- `reversalPotential`: reversal Nernst potential of the ion, Eion.
- `numberOfGatings`: number of gating variable(s) of the ionic current.
- `exponents`: exponent(s) of the gating variable(s), a1 (and a2).
- `activationSteadyStateGating`: equilibrium function of the activation gating variable.
- `activationTimeConstant`: time constant function/constant of the activation gating variable.
- `inactivationSteadyStateGating`: equilibrium function of the inactivation gating variable.
- `inactivationTimeConstant`: time constant function/constant of the inactivation gating variable.
- `calciumDependency`: flag indicating if the ionic current depends on calcium.
- `MgDependency`: flag indicating if the ionic current depends on magnesium.


# Example
```jldoctest
julia> NaCurrent = initializeCurrent("Na", VNa, numberOfGatings=2, exponents=[3, 1],
    activationSteadyStateGating=mNa_inf, activationTimeConstant=tau_mNa,
    inactivationSteadyStateGating=hNa_inf, inactivationTimeConstant=tau_hNa)
IonCurrent("Na", 2, [3, 1], Function[mNa_inf, hNa_inf], Function[tau_mNa, tau_hNa], 50.0, false, false)
```
"""
function initializeCurrent(name::String, reversalPotential::Union{Int64, Float64}; numberOfGatings::Int64=1, exponents::Union{Int64, Vector{Int64}}=1,
    activationSteadyStateGating::Union{Function, Float64}=0., activationTimeConstant::Union{Function, Float64, Int64}=NaN,
    inactivationSteadyStateGating::Union{Function, Float64}=0., inactivationTimeConstant::Union{Function, Float64, Int64}=NaN,
    calciumDependency::Bool=false, MgDependency::Bool=false)

    # If no activation steady state function is provided, throw error
    if !(typeof(activationSteadyStateGating) <: Function)
        error("Please enter a steady state activation gating value of type Function!")
    end

    # If no activation time constant is provided, throw error
    if !(typeof(activationTimeConstant) <: Function) && isnan(activationTimeConstant)
        error("Please enter an activation time constant of type Int64, Float64 or Function!")
    end

    # Build the activation time constant function
    activationTau_Function = transformToFunction(activationTimeConstant)

    # If we only have an activation gating variable
    if numberOfGatings == 1
        # If more than one exponent is provided, throw error
        if !isa(exponents, Int64)
            error("Given that you only have one gating variable in the ionic current, exponents must be of type Int64!")
        end

        # Create the IonCurrent structure
        return IonCurrent(name, numberOfGatings, [exponents, ], [activationSteadyStateGating, ],
            [activationTau_Function, ], Float64(reversalPotential), calciumDependency, MgDependency)

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

        # If no inactivation time constant is provided, throw error
        if !(typeof(inactivationTimeConstant) <: Function) && isnan(inactivationTimeConstant)
            error("Please enter an inactivation time constant of type Int64, Float64 or Function!")
        end

        # Build the activation time constant function
        inactivationTau_Function = transformToFunction(inactivationTimeConstant)

        # Create the IonCurrent structure
        return IonCurrent(name, numberOfGatings, exponents, [activationSteadyStateGating, inactivationSteadyStateGating],
            [activationTau_Function, inactivationTau_Function], Float64(reversalPotential), calciumDependency, MgDependency)

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
            globalMgDependency = true
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
        Float64(reversaleLeakagePotential), calciumDynamics, globalCalciumDependency, globalMgDependency)
end


