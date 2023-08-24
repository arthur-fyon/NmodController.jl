module NmodController
# Calling dependencies
using SymPy

# Calling all the files of the package
include("mathFunctions.jl")
include("neuron.jl")
include("elementsOfS.jl")
include("DIC.jl")

# Export meaningful functions
export boltz
export tauX
export IonCurrent
export CalciumDynamic
export NeuronCB
export initializeCurrent
export initializeCalciumDynamics
export initializeNeuronModel
export computeUnweightedUnigatedVoltageElementOfS
export computeUnweightedBigated1VoltageElementOfS
export computeUnweightedBigated2VoltageElementOfS
export computeUnweightedUnigatedCaElementOfS
export computeUnweightedBigated1CaElementOfS
export computeUnweightedBigated2CaElementOfS
export computeUnweightedUnigatedMgElementOfS
export computeUnweightedBigated1MgElementOfS
export computeEquilibriumOfCalcium
export computeEquilibriumCaCurrent
export computeWeightedElementOfS
export computeWeightedVoltageElementOfS
export computeWeightedCaElementOfS
export computeWeightedCaElementOfS
export computeWeights
export logDist
export computeS
export computeDICs

end
