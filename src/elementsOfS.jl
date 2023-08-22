# Functions that compute only voltage-dependent elements of the sensitivity matrix
# Only one gating variable
function computeUnweightedUnigatedVoltageElementOfS(ionCurrent::IonCurrent)
    # Compute the unweighted element of S
    exp1 = ionCurrent.exponents[1]
    derivateOfGating1 = lambdify(SymPy.diff(ionCurrent.steadyStateGatings[1]))
    Sij(V) = exp1*ionCurrent.steadyStateGatings[1](V)^(exp1-1) * (V-ionCurrent.reversalPotential) * derivateOfGating1(V)
    return Sij
end

# Two gating variable and derivatation according to the activation
function computeUnweightedBigated1VoltageElementOfS(ionCurrent::IonCurrent)
    # Compute the unweighted element of S
    exp1 = ionCurrent.exponents[1]
    exp2 = ionCurrent.exponents[2]
    derivateOfGating1 = lambdify(SymPy.diff(ionCurrent.steadyStateGatings[1]))
    Sij(V) = exp1*ionCurrent.steadyStateGatings[1](V)^(exp1-1) * ionCurrent.steadyStateGatings[2](V)^exp2 * (V-ionCurrent.reversalPotential) * derivateOfGating1(V)
    return Sij
end

# Two gating variable and derivatation according to the inactivation
function computeUnweightedBigated2VoltageElementOfS(ionCurrent::IonCurrent)
    # Compute the unweighted element of S
    exp1 = ionCurrent.exponents[1]
    exp2 = ionCurrent.exponents[2]
    derivateOfGating2 = lambdify(SymPy.diff(ionCurrent.steadyStateGatings[2]))
    Sij(V) = ionCurrent.steadyStateGatings[1](V)^exp1 * exp2*ionCurrent.steadyStateGatings[2](V)^(exp2-1) * (V-ionCurrent.reversalPotential) * derivateOfGating2(V)
    return Sij
end

# Function that compute Ca-dependent elements of the sensitivity matrix
# Only one gating variable
function computeUnweightedUnigatedCaElementOfS(ionCurrent::IonCurrent, CaEquilibrium::Function)
    # Compute the part of the unweighted element of S without derivating calcium
    exp1 = ionCurrent.exponents[1]
    partialDerivate1OfGating1 = lambdify(partialDiff1(ionCurrent.steadyStateGatings[1]))
    Sij1(V) = exp1*ionCurrent.steadyStateGatings[1](V, CaEquilibrium(V))^(exp1-1) * (V-ionCurrent.reversalPotential) * 
        partialDerivate1OfGating1(V, CaEquilibrium(V))

    # Compute the part of the unweighted element of S by derivating only voltage-dependent calcium
    exp1 = ionCurrent.exponents[1]
    partialDerivate2OfGating1 = lambdify(partialDiff2(ionCurrent.steadyStateGatings[1]))
    derivativeOfCa = lambdify(SymPy.diff(CaEquilibrium))
    Sij2(V) = exp1*ionCurrent.steadyStateGatings[1](V, CaEquilibrium(V))^(exp1-1) * (V-ionCurrent.reversalPotential) * 
        partialDerivate2OfGating1(V, CaEquilibrium(V)) * derivativeOfCa(V)

    return Sij1, Sij2
end

# Two gating variable and derivatation according to the activation
function computeUnweightedBigated1CaElementOfS(ionCurrent::IonCurrent, CaEquilibrium::Function)
    # Compute the part of the unweighted element of S without derivating calcium
    exp1 = ionCurrent.exponents[1]
    exp2 = ionCurrent.exponents[2]
    partialDerivate1OfGating1 = lambdify(partialDiff1(ionCurrent.steadyStateGatings[1]))
    Sij1(V) = exp1*ionCurrent.steadyStateGatings[1](V, CaEquilibrium(V))^(exp1-1) * 
        steadyStateGatings[2](V)^exp2 * (V-ionCurrent.reversalPotential) * 
        partialDerivate1OfGating1(V, CaEquilibrium(V))

    # Compute the part of the unweighted element of S by derivating only voltage-dependent calcium
    exp1 = ionCurrent.exponents[1]
    exp2 = ionCurrent.exponents[2]
    partialDerivate2OfGating1 = lambdify(partialDiff2(ionCurrent.steadyStateGatings[1]))
    derivativeOfCa = lambdify(SymPy.diff(CaEquilibrium))
    Sij2(V) = exp1*ionCurrent.steadyStateGatings[1](V, CaEquilibrium(V))^(exp1-1) *
        steadyStateGatings[2](V)^exp2 * (V-ionCurrent.reversalPotential) * 
        partialDerivate2OfGating1(V, CaEquilibrium(V)) * derivativeOfCa(V)

    return Sij1, Sij2
end

# Two gating variable and derivatation according to the inactivation
function computeUnweightedBigated2CaElementOfS(ionCurrent::IonCurrent, CaEquilibrium::Function)
    # Compute the part of the unweighted element of S without derivating calcium
    exp1 = ionCurrent.exponents[1]
    exp2 = ionCurrent.exponents[2]
    partialDerivate1OfGating2 = lambdify(partialDiff1(ionCurrent.steadyStateGatings[2]))
    Sij1(V) = ionCurrent.steadyStateGatings[1](V, CaEquilibrium(V))^exp1 * 
        exp2*steadyStateGatings[2](V)^(exp2-1) * (V-ionCurrent.reversalPotential) * 
        partialDerivate1OfGating2(V, CaEquilibrium(V))

    # Compute the part of the unweighted element of S by derivating only voltage-dependent calcium
    exp1 = ionCurrent.exponents[1]
    exp2 = ionCurrent.exponents[2]
    partialDerivate2OfGating2 = lambdify(partialDiff2(ionCurrent.steadyStateGatings[2]))
    derivativeOfCa = lambdify(SymPy.diff(CaEquilibrium))
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
    derivateOfGating1 = lambdify(partialDiff1(ionCurrent.steadyStateGatings[1]))
    Sij(V) = exp1*ionCurrent.steadyStateGatings[1](V)^(exp1-1) * (V-ionCurrent.reversalPotential) * derivateOfGating1(V, Mg)
    return Sij
end

# Two gating variable and derivatation according to the activation
function computeUnweightedBigated1MgElementOfS(ionCurrent::IonCurrent, Mg::Float64)
    # Compute the unweighted element of S
    exp1 = ionCurrent.exponents[1]
    exp2 = ionCurrent.exponents[2]
    derivateOfGating1 = lambdify(partialDiff1(ionCurrent.steadyStateGatings[1]))
    Sij(V) = exp1*ionCurrent.steadyStateGatings[1](V)^(exp1-1) * ionCurrent.steadyStateGatings[2](V)^exp2 * (V-ionCurrent.reversalPotential) * derivateOfGating1(V, Mg)
    return Sij
end

# Two gating variable and derivatation according to the inactivation
function computeUnweightedBigated2MgElementOfS(ionCurrent::IonCurrent, Mg::Float64)
    # Compute the unweighted element of S
    exp1 = ionCurrent.exponents[1]
    exp2 = ionCurrent.exponents[2]
    derivateOfGating2 = lambdify(partialDiff1(ionCurrent.steadyStateGatings[2]))
    Sij(V) = ionCurrent.steadyStateGatings[1](V)^exp1 * exp2*ionCurrent.steadyStateGatings[2](V)^(exp2-1) * (V-ionCurrent.reversalPotential) * derivateOfGating2(V, Mg)
    return Sij
end

# Function that compute the equilibrium value of the calcium
function computeEquilibriumOfCalcium(neuron::NeuronCB)
    # If no calcium dependency in the model, return 0
    if !neuron.globalCalciumDependency
        return 0.
    end

    # Extracting calcium dynamics
    calciumDynamics = neuron.calciumDynamics

    # Computing the value of calcium at equilibrium depending on the number of current in the calcium dynamics
    if length(calciumDynamics.coefficients) == 1
        # Extracting the current of interest
        currentIndex = findfirst(x -> occursin(calciumDynamics.currentNames[1], x.name), neuron.ionCurrents)
        if isa(currentIndex, Nothing)
            error("Current of the calcium dynamics not found in the model!")
        end
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
            currentIndex = findfirst(x -> occursin(currentName, x.name), neuron.ionCurrents)
            if isa(currentIndex, Nothing)
                error("Current of the calcium dynamics not found in the model!")
            end
            currentIndices[i] = currentIndex
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

# Function that computes weighted elements of the sensitivity matrix depending on the calcium/Mg dependency of the current
function computeWeightedElementOfS(ionCurrent::IonCurrent, tauFast::Function, tauSlow::Function, tauUltraslow::Function, 
    CaEquilibrium::Union{Float64, Function}, tauCa::Union{Float64, Function}, Mg::Float64)
    # If there is no calcium nor Mg dependency in the current
    if !(ionCurrent.calciumDependency) && !(ionCurrent.MgDependency)
        # Compute weighted elements of S
        Sf, Ss, Su = computeWeightedVoltageElementOfS(ionCurrent, tauFast, tauSlow, tauUltraslow)
        return Sf, Ss, Su
    # If there is only calcium dependency in the ion current
    elseif ionCurrent.calciumDependency && !(ionCurrent.MgDependency)
        # If no tauCa in argument, throw error
        if isnan(tauCa)
            error("Please enter a calcium time constant of type Float64 or Function!")
        end

        # Build the calcium time constant function
        tauCa_Function = transformToFunction(tauCa)

        # Compute weighted elements of S
        SfCa, SsCa, SuCa = computeWeightedCaElementOfS(ionCurrent, tauFast, tauSlow, tauUltraslow, CaEquilibrium, tauCa_Function)
        return SfCa, SsCa, SuCa
    # If there is only Mg dependency in the model
    elseif !(ionCurrent.calciumDependency) && ionCurrent.MgDependency
        # If no Mg in argument, throw error
        if isnan(Mg)
            error("Please enter a Mg value of type Float64!")
        end

        # Compute weighted elements of S
        SfMg, SsMg, SuMg = computeWeightedMgElementOfS(ionCurrent, tauFast, tauSlow, tauUltraslow, Mg)
        return SfMg, SsMg, SuMg
    else
        error("Ion current that depends both on calcium and Mg are not taken into account, sorry!")
    end
end

# Function that computes weighted elements of the sensitivity matrix that are only voltage dependent
function computeWeightedVoltageElementOfS(ionCurrent::IonCurrent, tauFast::Function, tauSlow::Function, tauUltraslow::Function)
    # If the current only has 1 gating variable
    if ionCurrent.numberOfGatings == 1
        # Compute the weights
        wfs, wsu = computeWeights(ionCurrent.timeConstants[1], tauFast, tauSlow, tauUltraslow)

        # Compute the unweighted element of S
        Sij = computeUnweightedUnigatedVoltageElementOfS(ionCurrent)

        # Compute the 3 weighted elements of S
        Sf1(V) = ionCurrent.steadyStateGatings[1](V)^ionCurrent.exponents[1] + wfs(V) * Sij(V)
        Ss1(V) = (wsu(V) - wfs(V)) * Sij(V)
        Su1(V) = (1 - wsu(V)) * Sij(V)
        return Sf1, Ss1, Su1
    # Otherwise the current has 2 gating variables
    else
        # Compute the weights
        wfs1, wsu1 = computeWeights(ionCurrent.timeConstants[1], tauFast, tauSlow, tauUltraslow)
        wfs2, wsu2 = computeWeights(ionCurrent.timeConstants[2], tauFast, tauSlow, tauUltraslow)

        # Compute the unweighted elements of S
        Sij1 = computeUnweightedBigated1VoltageElementOfS(ionCurrent)
        Sij2 = computeUnweightedBigated2VoltageElementOfS(ionCurrent)

        # Compute the 3 weighted elements of S
        Sf2(V) = ionCurrent.steadyStateGatings[1](V)^ionCurrent.exponents[1]*
            ionCurrent.steadyStateGatings[2](V)^ionCurrent.exponents[2] + wfs1(V) * Sij1(V) + wfs2(V) * Sij2(V)
        Ss2(V) = (wsu1(V) - wfs1(V)) * Sij1(V) + (wsu2(V) - wfs2(V)) * Sij2(V)
        Su2(V) = (1 - wsu1(V)) * Sij1(V) + (1 - wsu2(V)) * Sij2(V)
        return Sf2, Ss2, Su2
    end
end

# Function that computes weighted elements of the sensitivity matrix that are calcium dependent
function computeWeightedCaElementOfS(ionCurrent::IonCurrent, tauFast::Function, tauSlow::Function, tauUltraslow::Function, CaEquilibrium::Function, tauCa::Function)
    # If the current only has 1 gating variable
    if ionCurrent.numberOfGatings == 1
        # Compute the weights
        wfs, wsu = computeWeights(ionCurrent.timeConstants[1], tauFast, tauSlow, tauUltraslow)
        wfsCa, wsuCa = computeWeights(tauCa, tauFast, tauSlow, tauUltraslow)

        # Compute the unweighted elements of S
        Sij, SijCa = computeUnweightedUnigatedCaElementOfS(ionCurrent, CaEquilibrium)

        # Compute the 3 weighted elements of S
        Sf1(V) = ionCurrent.steadyStateGatings[1](V, CaEquilibrium(V))^ionCurrent.exponents[1] + wfs(V) * Sij(V) + wfsCa(V) * SijCa(V)
        Ss1(V) = (wsu(V) - wfs(V)) * Sij(V) + (wsuCa(V) - wfsCa(V)) * SijCa(V)
        Su1(V) = (1 - wsu(V)) * Sij(V) + (1 - wsuCa(V)) * SijCa(V)
        return Sf1, Ss1, Su1
    # Otherwise the current has 2 gating variables
    else
        # Compute the weights
        wfs1, wsu1 = computeWeights(ionCurrent.timeConstants[1], tauFast, tauSlow, tauUltraslow)
        wfs2, wsu2 = computeWeights(ionCurrent.timeConstants[2], tauFast, tauSlow, tauUltraslow)
        wfsCa, wsuCa = computeWeights(tauCa, tauFast, tauSlow, tauUltraslow)

        # Compute the unweighted elements of S
        Sij1, SijCa1 = computeUnweightedBigated1CaElementOfS(ionCurrent, CaEquilibrium)
        Sij2, SijCa2 = computeUnweightedBigated2CaElementOfS(ionCurrent, CaEquilibrium)

        # Compute the 3 weighted elements of S
        Sf2(V) = ionCurrent.steadyStateGatings[1](V, CaEquilibrium(V))^ionCurrent.exponents[1]*
            ionCurrent.steadyStateGatings[2](V, CaEquilibrium(V))^ionCurrent.exponents[2] + wfs1(V) * Sij1(V) + wfs2(V) * Sij2(V) + wfsCa(V) * SijCa1(V) + wfsCa(V) * SijCa2(V)
        Ss2(V) = (wsu1(V) - wfs1(V)) * Sij1(V) + (wsu2(V) - wfs2(V)) * Sij2(V) + (wsuCa(V) - wfsCa(V)) * SijCa1(V) + (wsuCa(V) - wfsCa(V)) * SijCa2(V)
        Su2(V) = (1 - wsu1(V)) * Sij1(V) + (1 - wsu2(V)) * Sij2(V) + (1 - wsuCa(V)) * SijCa1(V) + (1 - wsuCa(V)) * SijCa2(V)
        return Sf2, Ss2, Su2
    end
end

# Function that computes weighted elements of the sensitivity matrix that are Mg dependent
function computeWeightedCaElementOfS(ionCurrent::IonCurrent, tauFast::Function, tauSlow::Function, tauUltraslow::Function, Mg::Float64)
    # If the current only has 1 gating variable
    if ionCurrent.numberOfGatings == 1
        # Compute the weights
        wfs, wsu = computeWeights(ionCurrent.timeConstants[1], tauFast, tauSlow, tauUltraslow)

        # Compute the unweighted elements of S
        Sij = computeUnweightedUnigatedMgElementOfS(ionCurrent, Mg)

        # Compute the 3 weighted elements of S
        Sf1(V) = ionCurrent.steadyStateGatings[1](V, Mg)^ionCurrent.exponents[1] + wfs(V) * Sij(V)
        Ss1(V) = (wsu(V) - wfs(V)) * Sij(V)
        Su1(V) = (1 - wsu(V)) * Sij(V)
        return Sf1, Ss1, Su1
    # Otherwise the current has 2 gating variables
    else
        # Compute the weights
        wfs1, wsu1 = computeWeights(ionCurrent.timeConstants[1], tauFast, tauSlow, tauUltraslow)
        wfs2, wsu2 = computeWeights(ionCurrent.timeConstants[2], tauFast, tauSlow, tauUltraslow)

        # Compute the unweighted elements of S
        Sij1 = computeUnweightedBigated1MgElementOfS(ionCurrent, Mg)
        Sij2 = computeUnweightedBigated2MgElementOfS(ionCurrent, Mg)

        # Compute the 3 weighted elements of S
        Sf(V) = ionCurrent.steadyStateGatings[1](V, Mg)^ionCurrent.exponents[1]*
            ionCurrent.steadyStateGatings[2](V, Mg)^ionCurrent.exponents[2] + wfs1(V) * Sij1(V) + wfs2(V) * Sij2(V)
        Ss(V) = (wsu1(V) - wfs1(V)) * Sij1(V) + (wsu2(V) - wfs2(V)) * Sij2(V)
        Su(V) = (1 - wsu1(V)) * Sij1(V) + (1 - wsu2(V)) * Sij2(V)
        return Sf, Ss, Su
    end
end
