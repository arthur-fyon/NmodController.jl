# Prebuild gating functions
boltz(V, A, B) = 1 / (1 + exp((V+A)/B))
tauX(V, A, B, D, E) = A - B / (1+exp((V+D)/E))

# Data structure that defines a neuronal current
struct ionCurrent
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
struct calciumDynamic
    numberOfCurrents::Int64
    currentNames::Vector{String}
    coefficients::Vector{Float64}
    nernstEquilibriumValue::Float64
    timeConstant::Float64
end

# Data structure that defines the neuronal model
struct neuronCB
    C::Float64
    numberOfIonCurrents::Int64
    ionCurrents::Vector{ionCurrent}
    calciumDynamics::calciumDynamic
    globalCalciumDependency::Bool
    globalMgDependency::Bool
end

# Function that initialize a certain type of current
function initializeCurrent(name::String, reversalPotential::Union{Int64, Float64}; numberOfGatings::Int64=1, exponents::Union{Int64, Vector{Int64}}=1,
    activationSteadyStateGating::Union{Function, Float64}=0., activationTimeConstants::Union{Function, Float64, Int64}=0.,
    inactivationSteadyStateGating::Union{Function, Float64}=0., inactivationTimeConstants::Union{Function, Float64, Int64}=0.,
    calciumDependency::Bool=false, MgDependency::Bool=false)

    # If no activation steady state function is provided, throw error
    if activationSteadyStateGating ≠ Function
        error("Please enter a steady state activation gating value of type Function")
    end

    # If we only have an activation gating variable
    if numberOfGatings == 1
        # If more than one exponent is provided, throw error
        if typeof(exponents) ≠ Int64
            error("Given that you only have one gating variable in the ionic current, exponents must be of type Int64")
        end

        # Create the ionCurrent structure
        return ionCurrent(name, numberOfGatings, [exponents, ], [activationSteadyStateGating, ],
            [activationTimeConstants, ], Float64(reversalPotential), calciumDependency, MgDependency)

    # If we have both an activation and inactivation gating variables
    elseif numberOfGatings == 2
        # If anything else than 2 exponents are provided, throw error
        if typeof(exponents) ≠ Vector{Int64} && length(exponents) ≠ 2
            error("Given that you have two gating variables in the ionic current, exponents must be of type Vector{Int64}")
        end
        # If no inactivation steady state function is provided, throw error
        if inactivationSteadyStateGating ≠ Function
            error("Please enter a steady state inactivation gating value of type Function")
        end

        # Create the ionCurrent structure
        return ionCurrent(name, numberOfGatings, exponents, [activationSteadyStateGating, inactivationSteadyStateGating],
            [activationTimeConstants, inactivationTimeConstants], Float64(reversalPotential), calciumDependency, MgDependency)

    # If less than 1 or more than 2 gating variables are provided, throw error
    else
        error("Number of gating variables must be equal to 1 or 2!")
    end
end

# Function that initialize calcium dynamics
function initializeCalciumDynamics(currentNames::Union{String, Vector{String}}, coefficients::Union{Int64, Float64, Vector{Int64}, Vector{Float64}},
    nernstEquilibriumValue::Union{Int64, Float64}, timeConstant::Union{Int64, Float64})

    # If there is only one current in the calcium dynamics
    if (typeof(coefficients) == Int64 || typeof(coefficients) == Float64) && typeof(currentNames) == String
        # Create the calciumDynamic structure
        return calciumDynamic(1, [currentNames, ], [Float64(coefficients), ], Float64(nernstEquilibriumValue), Float64(timeConstant))

    # If there is more than one current in the calcium dynamics
    elseif (typeof(coefficients) == Vector{Int64} || typeof(coefficients) == Vector{Float64}) && typeof(currentNames) == Vector{String}
        # If length of currentNames and coefficients are the same
        if length(coefficients) == length(currentNames)
            # Create the calciumDynamic structure
            return calciumDynamic(length(coefficients), currentNames, Vector{Float64}(coefficients), Float64(nernstEquilibriumValue), Float64(timeConstant))
        # If length of currentNames and coefficients are not the same, throw error
        else
            error("Length of currentNames and coefficients must be the same!")
        end

    # Else, throw an error
    else
        error("currentNames and coefficients must be both of type Vector or of type String and Float64/Int64 respectively!")
    end
end
