# Include partial derivative functions
include("mathFunctions.jl")

# Function that computes the DICs or the sensitivity matrix
function computeDICs(neuron::NeuronCB)

end

# Function that compute only voltage-dependent elements of the sensitivity matrix
function computeUnweightedVoltageElementOfS(ionCurrent::IonCurrent, n::Int64)
    # If the current only have one gating variable
    if ionCurrent.numberOfGatings == 1
        Sij(V, Ca) = ionCurrent.exponents[1]*ionCurrent.steadyStateGatings #work in progress
    end
end