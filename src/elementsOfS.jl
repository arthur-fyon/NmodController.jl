# Include partial derivative functions
include("mathFunctions.jl")

# Functions that compute only voltage-dependent elements of the sensitivity matrix
# Only one gating variable
function computeUnweightedUnigatedVoltageElementOfS(ionCurrent::IonCurrent)
    # Compute the unweighted element of S
    exp1 = ionCurrent.exponents[1]
    derivateOfGating1 = SymPy.diff(ionCurrent.steadyStateGatings[1])
    Sij(V) = exp1*ionCurrent.steadyStateGatings[1](V)^(exp1-1) * (V-ionCurrent.reversalPotential) * derivateOfGating1(V)
    return Sij
end

# Two gating variable and derivatation according to the activation
function computeUnweightedBigated1VoltageElementOfS(ionCurrent::IonCurrent)
    # Compute the unweighted element of S
    exp1 = ionCurrent.exponents[1]
    exp2 = ionCurrent.exponents[2]
    derivateOfGating1 = SymPy.diff(ionCurrent.steadyStateGatings[1])
    Sij(V) = exp1*ionCurrent.steadyStateGatings[1](V)^(exp1-1) * ionCurrent.steadyStateGatings[2](V)^exp2 * (V-ionCurrent.reversalPotential) * derivateOfGating1(V)
    return Sij
end

# Two gating variable and derivatation according to the inactivation
function computeUnweightedBigated2VoltageElementOfS(ionCurrent::IonCurrent)
    # Compute the unweighted element of S
    exp1 = ionCurrent.exponents[1]
    exp2 = ionCurrent.exponents[2]
    derivateOfGating2 = SymPy.diff(ionCurrent.steadyStateGatings[2])
    Sij(V) = ionCurrent.steadyStateGatings[1](V)^exp1 * exp2*ionCurrent.steadyStateGatings[2](V)^(exp2-1) * (V-ionCurrent.reversalPotential) * derivateOfGating2(V)
    return Sij
end

# Function that compute Ca-dependent elements of the sensitivity matrix
# Only one gating variable
function computeUnweightedUnigatedCaElementOfS(ionCurrent::IonCurrent, calciumDynamics::CalciumDynamic, neuron::NeuronCB)
    # First computing the value of calcium at equilibrium
    CaEquilibrium = computeEquilibriumOfCalcium(calciumDynamics, neuron)

    # Compute the part of the unweighted element of S without derivating calcium
    exp1 = ionCurrent.exponents[1]
    partialDerivate1OfGating1 = partialDiff1(ionCurrent.steadyStateGatings[1])
    Sij1(V) = exp1*ionCurrent.steadyStateGatings[1](V, CaEquilibrium(V))^(exp1-1) * (V-ionCurrent.reversalPotential) * 
        partialDerivate1OfGating1(V, CaEquilibrium(V))

    # Compute the part of the unweighted element of S by derivating only voltage-dependent calcium
    exp1 = ionCurrent.exponents[1]
    partialDerivate2OfGating1 = partialDiff2(ionCurrent.steadyStateGatings[1])
    derivativeOfCa = SymPy.diff(CaEquilibrium)
    Sij2(V) = exp1*ionCurrent.steadyStateGatings[1](V, CaEquilibrium(V))^(exp1-1) * (V-ionCurrent.reversalPotential) * 
        partialDerivate2OfGating1(V, CaEquilibrium(V)) * derivativeOfCa(V)

    return Sij1, Sij2
end

# Two gating variable and derivatation according to the activation
function computeUnweightedBigated1CaElementOfS(ionCurrent::IonCurrent, calciumDynamics::CalciumDynamic, neuron::NeuronCB)
    # First computing the value of calcium at equilibrium
    CaEquilibrium = computeEquilibriumOfCalcium(calciumDynamics, neuron)

    # Compute the part of the unweighted element of S without derivating calcium
    exp1 = ionCurrent.exponents[1]
    exp2 = ionCurrent.exponents[2]
    partialDerivate1OfGating1 = partialDiff1(ionCurrent.steadyStateGatings[1])
    Sij1(V) = exp1*ionCurrent.steadyStateGatings[1](V, CaEquilibrium(V))^(exp1-1) * 
        steadyStateGatings[2](V)^exp2 * (V-ionCurrent.reversalPotential) * 
        partialDerivate1OfGating1(V, CaEquilibrium(V))

    # Compute the part of the unweighted element of S by derivating only voltage-dependent calcium
    exp1 = ionCurrent.exponents[1]
    exp2 = ionCurrent.exponents[2]
    partialDerivate2OfGating1 = partialDiff2(ionCurrent.steadyStateGatings[1])
    derivativeOfCa = SymPy.diff(CaEquilibrium)
    Sij2(V) = exp1*ionCurrent.steadyStateGatings[1](V, CaEquilibrium(V))^(exp1-1) *
        steadyStateGatings[2](V)^exp2 * (V-ionCurrent.reversalPotential) * 
        partialDerivate2OfGating1(V, CaEquilibrium(V)) * derivativeOfCa(V)

    return Sij1, Sij2
end

# Two gating variable and derivatation according to the inactivation
function computeUnweightedBigated2CaElementOfS(ionCurrent::IonCurrent, calciumDynamics::CalciumDynamic, neuron::NeuronCB)
    # First computing the value of calcium at equilibrium
    CaEquilibrium = computeEquilibriumOfCalcium(calciumDynamics, neuron)

    # Compute the part of the unweighted element of S without derivating calcium
    exp1 = ionCurrent.exponents[1]
    exp2 = ionCurrent.exponents[2]
    partialDerivate1OfGating2 = partialDiff1(ionCurrent.steadyStateGatings[2])
    Sij1(V) = ionCurrent.steadyStateGatings[1](V, CaEquilibrium(V))^exp1 * 
        exp2*steadyStateGatings[2](V)^(exp2-1) * (V-ionCurrent.reversalPotential) * 
        partialDerivate1OfGating2(V, CaEquilibrium(V))

    # Compute the part of the unweighted element of S by derivating only voltage-dependent calcium
    exp1 = ionCurrent.exponents[1]
    exp2 = ionCurrent.exponents[2]
    partialDerivate2OfGating2 = partialDiff2(ionCurrent.steadyStateGatings[2])
    derivativeOfCa = SymPy.diff(CaEquilibrium)
    Sij2(V) = ionCurrent.steadyStateGatings[1](V, CaEquilibrium(V))^exp1 *
        exp2*steadyStateGatings[2](V)^(exp2-1) * (V-ionCurrent.reversalPotential) * 
        partialDerivate2OfGating2(V, CaEquilibrium(V)) * derivativeOfCa(V)

    return Sij1, Sij2
end

# Function that compute Mg-dependent elements of the sensitivity matrix
# Only one gating variable
function computeUnweightedUnigatedMgElementOfS(ionCurrent::IonCurrent, Mg::Float64)
    # Compute the unweighted element of S
    exp1 = ionCurrent.exponents[1]
    derivateOfGating1 = partialDiff1(ionCurrent.steadyStateGatings[1])
    Sij(V) = exp1*ionCurrent.steadyStateGatings[1](V)^(exp1-1) * (V-ionCurrent.reversalPotential) * derivateOfGating1(V, Mg)
    return Sij
end

# Two gating variable and derivatation according to the activation
function computeUnweightedBigated1MgElementOfS(ionCurrent::IonCurrent, Mg::Float64)
    # Compute the unweighted element of S
    exp1 = ionCurrent.exponents[1]
    exp2 = ionCurrent.exponents[2]
    derivateOfGating1 = partialDiff1(ionCurrent.steadyStateGatings[1])
    Sij(V) = exp1*ionCurrent.steadyStateGatings[1](V)^(exp1-1) * ionCurrent.steadyStateGatings[2](V)^exp2 * (V-ionCurrent.reversalPotential) * derivateOfGating1(V, Mg)
    return Sij
end

# Two gating variable and derivatation according to the inactivation
function computeUnweightedBigated1MgElementOfS(ionCurrent::IonCurrent, Mg::Float64)
    # Compute the unweighted element of S
    exp1 = ionCurrent.exponents[1]
    exp2 = ionCurrent.exponents[2]
    derivateOfGating2 = partialDiff1(ionCurrent.steadyStateGatings[2])
    Sij(V) = ionCurrent.steadyStateGatings[1](V)^exp1 * exp2*ionCurrent.steadyStateGatings[2](V)^(exp2-1) * (V-ionCurrent.reversalPotential) * derivateOfGating2(V)
    return Sij
end

# Function that compute the equilibrium value of the calcium
function computeEquilibriumOfCalcium(calciumDynamics::CalciumDynamic, neuron::NeuronCB)
    # Computing the value of calcium at equilibrium depending on the number of current in the calcium dynamics
    if length(calciumDynamics.coefficients) == 1
        # Extracting the current of interest
        currentIndex = findfirst(x -> occursin(calciumDynamics.currentNames[1], x.name), neuron.ionCurrents)
        ionCurrentCa = neuron.ionCurrents[currentIndex]

        # Extracting the value of the maximal conductance and set it to 10 if not defined
        g = neuron.maximumConductances[currentIndex]
        if isnan(g)
            g = 10.
        end

        # Computing the equilibrium calcium current value
        currentValue = computeEquilibriumCaCurrent(calciumDynamics.coefficients[1], ionCurrentCa, g)

        # Compute the calcium equilibrium value
        CaEquilibrium1(V) = calciumDynamics.nernstEquilibriumValue + currentValue(V)
        return CaEquilibrium1

    elseif length(calciumDynamics.coefficients) == 2
        # Extracting the currents of interest
        currentIndices = zeros(Int64, 2)
        for (i, currentName) in enumerate(calciumDynamics.currentNames)
            currentIndices[i] = findfirst(x -> occursin(currentName, x.name), neuron.ionCurrents)
        end
        ionCurrentCa1 = neuron.ionCurrents[currentIndices[1]]
        ionCurrentCa2 = neuron.ionCurrents[currentIndices[2]]

        # Extracting the values of the maximal conductances and set it to 10 if not defined
        g1 = neuron.maximumConductances[currentIndices[1]]
        if isnan(g1)
            g1 = 10.
        end

        g2 = neuron.maximumConductances[currentIndices[2]]
        if isnan(g2)
            g2 = 10.
        end

        # Computing the equilibrium calcium current values
        currentValue1 = computeEquilibriumCaCurrent(calciumDynamics.coefficients[1], ionCurrentCa1, g1)
        currentValue2 = computeEquilibriumCaCurrent(calciumDynamics.coefficients[2], ionCurrentCa2, g2)

        CaEquilibrium2(V) = calciumDynamics.nernstEquilibriumValue + currentValue1(V) + currentValue2(V)
        return CaEquilibrium2
    else
        error("Calcium dynamics with more than 2 currents are not yet taken into account, sorry!")
    end
end

# Function that compute the equilibrium value of calcium current
function computeEquilibriumCaCurrent(coefficient::Float64, ionCurrent::IonCurrent, g::Float64)
    # Compute the current equilibrium values depending on the number of gating variables
    if ionCurrent.numberOfGatings == 1
        currentValue1(V) = coefficient * g * ionCurrent.steadyStateGatings[1](V)^ionCurrent.exponents[1] * (V-ionCurrent.reversalPotential)
        return currentValue1
    else
        currentValue2(V) = coefficient * g * ionCurrent.steadyStateGatings[1](V)^ionCurrent.exponents[1] * 
            ionCurrent.steadyStateGatings[2](V)^ionCurrent.exponents[2] * (V-ionCurrent.reversalPotential)
        return currentValue2
    end
end