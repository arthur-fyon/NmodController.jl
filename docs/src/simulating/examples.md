# Example on existing models in the litterature

## Simulating the uncontrolled STG model
Once the STG model is initialized, simulated it can be simplified using *NmodController.jl*. At first, we will simulate the model `STG_spiking` that exhibits a tonic spiking behavior with
- sodium current: $\bar{g}_\mathrm{Na} =$ 4000;
- T-type calcium current: $\bar{g}_\mathrm{CaT} =$ 3;
- slow calcium current: $\bar{g}_\mathrm{CaS} =$ 4;
- A-type potassium current: $\bar{g}_\mathrm{A} =$ 175;
- calcium controlled potassium current: $\bar{g}_\mathrm{KCa} =$ 110;
- delayed rectified potassium current: $\bar{g}_\mathrm{Kd} =$ 137;
- H type current: $\bar{g}_\mathrm{H} =$ 0.3;
- leakage current: $g_\mathrm{leak} =$ 0.01.

The next few lines of code show how to simulate the spiking STG model using *NmodController.jl*.

```julia
# First writing uncontrolled ODE function file
writeUncontrolledODEs(STG_spiking, filename="STG_ODE.jl")

# Including the newly written file
include("STG_ODE.jl")

# Definition of simulation time (in ms)
Tfinal = 2000
tspan  = (0., Tfinal)

# Definition of membrane capacitance and maximal conductance values
gNa   = 4000. # Sodium current maximal conductance
gCaT  = 3. # T-type calcium current maximal conductance
gCaS  = 4. # Slow calcium current maximal conductance
gA    = 175. # A-type potassium current maximal conductance
gKCa  = 110. # Calcium-activated potassium current maximal conductance
gKd   = 137. # Delayed-rectifier potassium current maximal conductance
gH    = 0.3 # H-current maximal conductance
gleak = 0.01 # Leak current maximal conductance
C = 1. # Membrane capacitance

# Input current definition (always a function of time)
Iapp(t) = 0.

# Parameter vector for simulations
p = (Iapp, gNa, gCaT, gCaS, gA, gKCa, gKd, gH, gleak, C)

# Initial conditions
V0  = -70.
Ca0 = 0.5
x0  = [V0, STG_mNa_inf(V0), STG_hNa_inf(V0), STG_mCaT_inf(V0), STG_hCaT_inf(V0), STG_mCaS_inf(V0), 
      STG_hCaS_inf(V0), STG_mA_inf(V0), STG_hA_inf(V0), STG_mKCa_inf(V0, Ca0), STG_mKd_inf(V0), STG_mH_inf(V0), Ca0]

# Simulation
prob = ODEProblem(STG_ODE, x0, tspan, p) # Describing the problem
sol  = solve(prob) # Solving the problem

# Retrieving variables
tt        = 0. : 0.2 : Tfinal
x         = sol(tt)
V_plot   = x[1, :]

# Plot
p = plot(tt/1e3, V_plot, xlims=(0, 2), xticks=[0, 2], ylims=(-100, 60), yticks=[-100, 60], linewidth=1.5, 
          legend=false, size=(600, 200), color=:gray30, margins=20px, dpi=500)
ylabel!("V (mV)")
xlabel!("t (s)")
```
![](https://raw.githubusercontent.com/arthur-fyon/NmodController.jl/main/docs/src/assets/STG_simu_spiking.png)

with `STG_ODE.jl` containing

```julia
#=
	This file contains differential equations describing the CB model of interest
=#

# Function that outputs values of variables derivatives
function STG_ODE(dx, x, p, t)
	# Parameters
	Iapp = p[1](t) # Time dependent applied current
	gNa = p[2] # Maximum conductance of current Na
	gCaT = p[3] # Maximum conductance of current CaT
	gCaS = p[4] # Maximum conductance of current CaS
	gA = p[5] # Maximum conductance of current A
	gKCa = p[6] # Maximum conductance of current KCa
	gKd = p[7] # Maximum conductance of current Kd
	gH = p[8] # Maximum conductance of current H
	gleak = p[9] # Leakage conductance
	C = p[10] # Membrane capacitance

	# Variables
	V = x[1] # Membrane voltage
	mNa = x[2] # Activation gating variable of current Na
	hNa = x[3] # Inactivation gating variable of current Na
	mCaT = x[4] # Activation gating variable of current CaT
	hCaT = x[5] # Inactivation gating variable of current CaT
	mCaS = x[6] # Activation gating variable of current CaS
	hCaS = x[7] # Inactivation gating variable of current CaS
	mA = x[8] # Activation gating variable of current A
	hA = x[9] # Inactivation gating variable of current A
	mKCa = x[10] # Activation gating variable of current KCa
	mKd = x[11] # Activation gating variable of current Kd
	mH = x[12] # Activation gating variable of current H
	Ca = x[13] # Ca concentration

	# ODEs
	dx[1] = (1/C) * (- gNa*mNa^3*hNa^1*(V - 50.0) -
					   gCaT*mCaT^3*hCaT^1*(V - 80.0) -
					   gCaS*mCaS^3*hCaS^1*(V - 80.0) -
					   gA*mA^3*hA^1*(V - -80.0) -
					   gKCa*mKCa^4*(V - -80.0) -
					   gKd*mKd^4*(V - -80.0) -
					   gH*mH^1*(V - -20.0) -
					   gleak*(V - -50.0) + Iapp)
	dx[2] = (1/STG_tau_mNa(V)) * (STG_mNa_inf(V) - mNa)
	dx[3] = (1/STG_tau_hNa(V)) * (STG_hNa_inf(V) - hNa)
	dx[4] = (1/STG_tau_mCaT(V)) * (STG_mCaT_inf(V) - mCaT)
	dx[5] = (1/STG_tau_hCaT(V)) * (STG_hCaT_inf(V) - hCaT)
	dx[6] = (1/STG_tau_mCaS(V)) * (STG_mCaS_inf(V) - mCaS)
	dx[7] = (1/STG_tau_hCaS(V)) * (STG_hCaS_inf(V) - hCaS)
	dx[8] = (1/STG_tau_mA(V)) * (STG_mA_inf(V) - mA)
	dx[9] = (1/STG_tau_hA(V)) * (STG_hA_inf(V) - hA)
	dx[10] = (1/STG_tau_mKCa(V)) * (STG_mKCa_inf(V, Ca) - mKCa)
	dx[11] = (1/STG_tau_mKd(V)) * (STG_mKd_inf(V) - mKd)
	dx[12] = (1/STG_tau_mH(V)) * (STG_mH_inf(V) - mH)
	dx[13] = (-0.94*gCaT*mCaT^3*hCaT^1*(V - 80.0) + -0.94*gCaS*mCaS^3*hCaS^1*(V - 80.0) + -Ca + 0.05) / 20.0
end
```

Note that a small adapation of `STG_ODE.jl` had been realized so that the name of the gating functions match the actual ones.

## Simulating the controlled STG model
As in Fyon et al., 2023 "Reliable neuromodulation from adaptive control of ion channel expression", the STG model slow and ultraslow DICs can be controlled by the slow calcium and A-type potassium currents to achieve a robust tonic spiking to bursting transition. 

The next few lines of code show how to simulate such controlled STG model using *NmodController.jl*.

```julia
# First writing controlled ODE function file
writeControlledODEs(STG_spiking, ["CaS", "A"], ["s", "u"], filename="ControlledSTG_ODE.jl")

# Including the newly written file
include("ControlledSTG_ODE.jl")

# Definition of simulation time (in ms)
Tfinal = 5000
tspan  = (0., Tfinal)

# Definition of membrane capacitance and maximal conductance values
gNa   = 4000. # Sodium current maximal conductance
gCaT  = 3. # T-type calcium current maximal conductance
gCaS  = 4. # Slow calcium current maximal conductance
gA    = 175. # A-type potassium current maximal conductance
gKCa  = 110. # Calcium-activated potassium current maximal conductance
gKd   = 137. # Delayed-rectifier potassium current maximal conductance
gH    = 0.3 # H-current maximal conductance
gleak = 0.01 # Leak current maximal conductance
C = 1. # Membrane capacitance

# Input current definition (always a function of time)
Iapp(t) = 0.

# Definition of controller parameters
α = 5e-3 # Rate of transfer between intracellular and membrane
β = 5e-3 # Rate of degradation of intracellular proteins
Kp = 3e-4 # Proprtional gain
Ki = 5e-6 # Integral gain
gsth(t) = 5. - 13 * (t > 500) # Reference gs(Vth)
guth(t) = 4. # Reference gu(Vth)

# Parameter vector for simulations
p = (Iapp, gNa, gCaT, gKCa, gKd, gH, gleak, C, α, β, Kp, Ki, STG_S_spiking(STG_Vth), gsth, guth)

# Initial conditions
V0  = -70.
Ca0 = 0.5
x0  = [V0, STG_mNa_inf(V0), STG_hNa_inf(V0), STG_mCaT_inf(V0), STG_hCaT_inf(V0), STG_mCaS_inf(V0), 
      STG_hCaS_inf(V0), STG_mA_inf(V0), STG_hA_inf(V0), STG_mKCa_inf(V0, Ca0), STG_mKd_inf(V0), STG_mH_inf(V0), Ca0,
      gCaS, gCaS, 0, gA, gA, 0]

# Simulation
prob = ODEProblem(ControlledSTG_ODE, x0, tspan, p) # Describing the problem
sol  = solve(prob) # Solving the problem

# Retrieving variables
tt        = 0. : 0.2 : Tfinal
x         = sol(tt)
V_plot    = x[1, :]
gCaS_plot = x[15, :]
gA_plot   = x[18, :]

# Plot
p1 = plot(tt/1e3, V_plot, xlims=(0, 5), xticks=[0, 5], ylims=(-100, 60), yticks=[-100, 60], linewidth=1.5, 
          legend=false, size=(600, 200), color=:gray30)
ylabel!("V (mV)")
p2 = plot(tt/1e3, gCaS_plot, xlims=(0, 5), xticks=[0, 5], ylims=(0, 20), yticks=[0, 20], linewidth=1.5, 
          legend=false, size=(600, 200), color=:gray30)
ylabel!("gCaS")
p3 = plot(tt/1e3, gA_plot, xlims=(0, 5), xticks=[0, 5], ylims=(0, 200), yticks=[0, 200], linewidth=1.5, 
          legend=false, size=(600, 200), color=:gray30)
ylabel!("gA")
xlabel!("t (s)")
CC = plot(p1, p2, p3, layout=(3, 1), size=(900, 600), margins=10px, dpi=500)
```
![](https://raw.githubusercontent.com/arthur-fyon/NmodController.jl/main/docs/src/assets/STG_simu_spiking_controlled.png)

with `ControlledSTG_ODE.jl` containing

```julia
#=
	This file contains differential equations describing the controlled CB model of interest
=#

# Function that outputs values of variables derivatives
function ControlledSTG_ODE(dx, x, p, t)
	# Parameters
	Iapp = p[1](t) # Time dependent applied current
	gNa = p[2] # Maximum conductance of current Na
	gCaT = p[3] # Maximum conductance of current CaT
	gKCa = p[4] # Maximum conductance of current KCa
	gKd = p[5] # Maximum conductance of current Kd
	gH = p[6] # Maximum conductance of current H
	gleak = p[7] # Leakage conductance
	C = p[8] # Membrane capacitance
	α = p[9] # Rate of transfer between intracellular and membrane
	β = p[10] # Rate of degradation of intracellular proteins
	Kp = p[11] # Proportional gain
	Ki = p[12] # Integral gain
	SVth = p[13] # Sensitivity matrix at threshold voltage
	gsth = p[14](t) # Reference gs(Vth)
	guth = p[15](t) # Reference gu(Vth)

	# Variables
	V = x[1] # Membrane voltage
	mNa = x[2] # Activation gating variable of current Na
	hNa = x[3] # Inactivation gating variable of current Na
	mCaT = x[4] # Activation gating variable of current CaT
	hCaT = x[5] # Inactivation gating variable of current CaT
	mCaS = x[6] # Activation gating variable of current CaS
	hCaS = x[7] # Inactivation gating variable of current CaS
	mA = x[8] # Activation gating variable of current A
	hA = x[9] # Inactivation gating variable of current A
	mKCa = x[10] # Activation gating variable of current KCa
	mKd = x[11] # Activation gating variable of current Kd
	mH = x[12] # Activation gating variable of current H
	Ca = x[13] # Ca concentration
	gCaSi = x[14] # Intracellular maximum conductance of current CaS
	gCaS = x[15] # Maximum conductance of current CaS
	zCaS = x[16] # Integral variable of current CaS
	gAi = x[17] # Intracellular maximum conductance of current A
	gA = x[18] # Maximum conductance of current A
	zA = x[19] # Integral variable of current A

	# ODEs
	dx[1] = (1/C) * (- gNa*mNa^3*hNa^1*(V - 50.0) -
					   gCaT*mCaT^3*hCaT^1*(V - 80.0) -
					   gCaS*mCaS^3*hCaS^1*(V - 80.0) -
					   gA*mA^3*hA^1*(V - -80.0) -
					   gKCa*mKCa^4*(V - -80.0) -
					   gKd*mKd^4*(V - -80.0) -
					   gH*mH^1*(V - -20.0) -
					   gleak*(V - -50.0) + Iapp)
	dx[2] = (1/STG_tau_mNa(V)) * (STG_mNa_inf(V) - mNa)
	dx[3] = (1/STG_tau_hNa(V)) * (STG_hNa_inf(V) - hNa)
	dx[4] = (1/STG_tau_mCaT(V)) * (STG_mCaT_inf(V) - mCaT)
	dx[5] = (1/STG_tau_hCaT(V)) * (STG_hCaT_inf(V) - hCaT)
	dx[6] = (1/STG_tau_mCaS(V)) * (STG_mCaS_inf(V) - mCaS)
	dx[7] = (1/STG_tau_hCaS(V)) * (STG_hCaS_inf(V) - hCaS)
	dx[8] = (1/STG_tau_mA(V)) * (STG_mA_inf(V) - mA)
	dx[9] = (1/STG_tau_hA(V)) * (STG_hA_inf(V) - hA)
	dx[10] = (1/STG_tau_mKCa(V)) * (STG_mKCa_inf(V, Ca) - mKCa)
	dx[11] = (1/STG_tau_mKd(V)) * (STG_mKd_inf(V) - mKd)
	dx[12] = (1/STG_tau_mH(V)) * (STG_mH_inf(V) - mH)
	dx[13] = (-0.94*gCaT*mCaT^3*hCaT^1*(V - 80.0) + -0.94*gCaS*mCaS^3*hCaS^1*(V - 80.0) + -Ca + 0.05) / 20.0

	# Retrieving which line of the sensitivity matrix matter
	timescales = [2, 3]

	# Retrieving which column of the sensitivity matrix belong to unmodulated conductances
	unmodulated = [1, 2, 5, 6, 7]

	# Retrieving which column of the sensitivity matrix belong to modulated conductances
	modulated = [3, 4]

	# Computing the right hand side of the linear system
	gDICr = [gsth, guth]
	gDICr = gDICr - SVth[timescales, unmodulated] * collect(p[2:6])

	# Computing the left hand side of the linear system
	Smod = SVth[timescales, modulated]

	# Computing the solution of the linear system
	g_r = \(Smod, gDICr)

	# Error signals and control inputs
	eCaS = g_r[1] - gCaS
	uCaS = Kp * eCaS + Ki * zCaS
	eA = g_r[2] - gA
	uA = Kp * eA + Ki * zA

	# ODEs of the controller
	dx[14] = α * gCaS - α * gCaSi - β * gCaSi + uCaS
	dx[15] = α * gCaSi - α * gCaS
	dx[16] = eCaS
	dx[17] = α * gA - α * gAi - β * gAi + uA
	dx[18] = α * gAi - α * gA
	dx[19] = eA
end
```

Note that a small adapation of `ControlledSTG_ODE.jl` had been realized so that the name of the gating functions match the actual ones.


## Simulating the uncontrolled DA model
Once the DA model is initialized, simulated it can be simplified using *NmodController.jl*. At first, we will simulate the model `DA_spiking` that exhibits a tonic spiking behavior with
- sodium current: $\bar{g}_\mathrm{Na} =$ 30;
- delayed rectified potassium current: $\bar{g}_\mathrm{Kd} =$ 5;
- L-type calcium current: $\bar{g}_\mathrm{CaL} =$ 0.03;
- N-type calcium current: $\bar{g}_\mathrm{CaN} =$ 0.03;
- ERG current: $\bar{g}_\mathrm{ERG} =$ 0.12;
- NMDA current: $\bar{g}_\mathrm{NMDA} =$ 0.1;
- leakage current: $g_\mathrm{leak} =$ 0.01.

The next few lines of code show how to simulate the spiking DA model using *NmodController.jl*.

```julia
# First writing uncontrolled ODE function file
writeUncontrolledODEs(DA_spiking, filename="DA_ODE.jl")

# Including the newly written file
include("DA_ODE.jl")

# Definition of simulation time (in ms)
Tfinal = 5000
tspan  = (0., Tfinal)

# Definition of membrane capacitance and maximal conductance values
gNa   = 30. # Sodium current maximal conductance
gKd   = 5. # Delayed-rectifier potassium current maximal conductance
gCaL  = 0.03 # L-type calcium current maximal conductance
gCaN  = 0.03 # N-type calcium current maximal conductance
gERG  = 0.12 # ERG current maximal conductance
gNMDA  = 0.1 # NMDA current maximal conductance
gleak = 0.01 # Leak current maximal conductance
C = 1. # Membrane capacitance
Mg = 1.4 # Magnesium concentration

# Input current definition (always a function of time)
Iapp(t) = 0.

# Parameter vector for simulations
p = (Iapp, gNa, gKd, gCaL, gCaN, gERG, gNMDA, gleak, C, Mg)

# Initial conditions
V0 = -90.
x0  = [V0, DA_mNa_inf(V0), DA_hNa_inf(V0), DA_mKd_inf(V0), DA_mCaL_inf(V0), DA_mCaN_inf(V0), 0., 0.]

# Simulation
prob = ODEProblem(DA_ODE, x0, tspan, p) # Describing the problem
sol  = solve(prob) # Solving the problem

# Retrieving variables
tt        = 0. : 0.2 : Tfinal
x         = sol(tt)
V_plot   = x[1, :]

# Plot
p = plot(tt/1e3, V_plot, xlims=(0, 5), xticks=[0, 5], ylims=(-100, 60), yticks=[-100, 60], linewidth=1.5, 
          legend=false, size=(600, 200), color=:gray30, margins=20px, dpi=500)
ylabel!("V (mV)")
xlabel!("t (s)")
```
![](https://raw.githubusercontent.com/arthur-fyon/NmodController.jl/main/docs/src/assets/DA_simu_spiking.png)

with `DA_ODE.jl` containing

```julia
#=
	This file contains differential equations describing the CB model of interest
=#

# Function that outputs values of variables derivatives
function DA_ODE(dx, x, p, t)
	# Parameters
	Iapp = p[1](t) # Time dependent applied current
	gNa = p[2] # Maximum conductance of current Na
	gKd = p[3] # Maximum conductance of current Kd
	gCaL = p[4] # Maximum conductance of current CaL
	gCaN = p[5] # Maximum conductance of current CaN
	gERG = p[6] # Maximum conductance of current ERG
	gNMDA = p[7] # Maximum conductance of current NMDA
	gleak = p[8] # Leakage conductance
	C = p[9] # Membrane capacitance
	Mg = p[10] # Mg concentration

	# Variables
	V = x[1] # Membrane voltage
	mNa = x[2] # Activation gating variable of current Na
	hNa = x[3] # Inactivation gating variable of current Na
	mKd = x[4] # Activation gating variable of current Kd
	mCaL = x[5] # Activation gating variable of current CaL
	mCaN = x[6] # Activation gating variable of current CaN
	oERG = x[7] # ERG potassium current activation
	iERG = x[8] # ERG current intermediate

	# ODEs
	dx[1] = (1/C) * (- gNa*mNa^3*hNa^1*(V - 60.0) -
					   gKd*mKd^3*(V - -85.0) -
					   gCaL*mCaL^2*(V - 60.0) -
					   gCaN*mCaN^1*(V - 60.0) -
					   gERG*oERG^1*(V - -85.0) -
					   gNMDA*DA_NMDA_inf(V, Mg)^1*(V - 0.0) -
					   gleak*(V - -50.0) + Iapp)
	dx[2] = (1/DA_tau_mNa(V)) * (DA_mNa_inf(V) - mNa)
	dx[3] = (1/DA_tau_hNa(V)) * (DA_hNa_inf(V) - hNa)
	dx[4] = (1/DA_tau_mKd(V)) * (DA_mKd_inf(V) - mKd)
	dx[5] = (1/DA_tau_mCaL(V)) * (DA_mCaL_inf(V) - mCaL)
	dx[6] = (1/DA_tau_mCaN(V)) * (DA_mCaN_inf(V) - mCaN)
	dx[7] = a0ERG(V) * (1-oERG-iERG) + biERG(V) * iERG - oERG * (aiERG(V)+b0ERG(V))
    dx[8] = aiERG(V) * oERG - biERG(V) * iERG
end

```

Note that a small adapation of `DA_ODE.jl` had been realized so that the name of the gating functions match the actual ones. Moreover, NMDA current had been changed to be instantaneaous (no gating variable dynamic) and the ERG gating variable dynamics were corrected since this current dynamic is non conventional.

## Simulating the controlled DA model
The DA model slow and ultraslow DICs can be controlled by the L-type and N-type calcium currents to achieve a robust tonic spiking to bursting transition. 

The next few lines of code show how to simulate such controlled DA model using *NmodController.jl*.

```julia
# First writing controlled ODE function file
writeControlledODEs(DA_spiking, ["CaL", "CaN"], ["s", "u"], filename="ControlledDA_ODE.jl")

# Including the newly written file
include("ControlledDA_ODE.jl")

# Definition of simulation time (in ms)
Tfinal = 10000
tspan  = (0., Tfinal)

# Definition of membrane capacitance and maximal conductance values
gNa   = 30. # Sodium current maximal conductance
gKd   = 5. # Delayed-rectifier potassium current maximal conductance
gCaL  = 0.03 # L-type calcium current maximal conductance
gCaN  = 0.03 # N-type calcium current maximal conductance
gERG  = 0.12 # ERG current maximal conductance
gNMDA  = 0.1 # NMDA current maximal conductance
gleak = 0.01 # Leak current maximal conductance
C = 1. # Membrane capacitance

# Input current definition (always a function of time)
Iapp(t) = 0.

# Definition of controller parameters
α = 5e-3 # Rate of transfer between intracellular and membrane
β = 5e-3 # Rate of degradation of intracellular proteins
Kp = 3e-4 # Proprtional gain
Ki = 5e-6 # Integral gain
gsth(t) = 0.5 - 4.5 * (t > 2000) # Reference gs(Vth)
guth(t) = 5. # Reference gu(Vth)

# Parameter vector for simulations
p = (Iapp, gNa, gKd, gERG, gNMDA, gleak, C, Mg, α, β, Kp, Ki, DA_S_spiking(DA_Vth), gsth, guth)

# Initial conditions
V0 = -90.
x0  = [V0, DA_mNa_inf(V0), DA_hNa_inf(V0), DA_mKd_inf(V0), DA_mCaL_inf(V0), DA_mCaN_inf(V0), 0., 0.,
      gCaL, gCaL, 0, gCaN, gCaN, 0]

# Simulation
prob = ODEProblem(ControlledDA_ODE, x0, tspan, p) # Describing the problem
sol  = solve(prob) # Solving the problem

# Retrieving variables
tt        = 0. : 0.2 : Tfinal
x         = sol(tt)
V_plot    = x[1, :]
gCaL_plot = x[10, :]
gCaN_plot   = x[13, :]

# Plot
p1 = plot(tt/1e3, V_plot, xlims=(0, 10), xticks=[0, 10], ylims=(-100, 60), yticks=[-100, 60], linewidth=1.5, 
          legend=false, size=(600, 200), color=:gray30)
ylabel!("V (mV)")
p2 = plot(tt/1e3, gCaL_plot, xlims=(0, 10), xticks=[0, 10], ylims=(0, 0.06), yticks=[0, 0.06], linewidth=1.5, 
          legend=false, size=(600, 200), color=:gray30)
ylabel!("gCaL")
p3 = plot(tt/1e3, gCaN_plot, xlims=(0, 10), xticks=[0, 10], ylims=(0, 0.2), yticks=[0, 0.2], linewidth=1.5, 
          legend=false, size=(600, 200), color=:gray30)
ylabel!("gCaN")
xlabel!("t (s)")
CC = plot(p1, p2, p3, layout=(3, 1), size=(900, 600), margins=10px, dpi=500)
```
![](https://raw.githubusercontent.com/arthur-fyon/NmodController.jl/main/docs/src/assets/DA_simu_spiking_controlled.png)

with `ControlledDA_ODE.jl` containing

```julia
#=
	This file contains differential equations describing the controlled CB model of interest
=#

# Function that outputs values of variables derivatives
function ControlledDA_ODE(dx, x, p, t)
	# Parameters
	Iapp = p[1](t) # Time dependent applied current
	gNa = p[2] # Maximum conductance of current Na
	gKd = p[3] # Maximum conductance of current Kd
	gERG = p[4] # Maximum conductance of current ERG
	gNMDA = p[5] # Maximum conductance of current NMDA
	gleak = p[6] # Leakage conductance
	C = p[7] # Membrane capacitance
	Mg = p[8] # Mg concentration
	α = p[9] # Rate of transfer between intracellular and membrane
	β = p[10] # Rate of degradation of intracellular proteins
	Kp = p[11] # Proportional gain
	Ki = p[12] # Integral gain
	SVth = p[13] # Sensitivity matrix at threshold voltage
	gsth = p[14](t) # Reference gs(Vth)
	guth = p[15](t) # Reference gu(Vth)

	# Variables
	V = x[1] # Membrane voltage
	mNa = x[2] # Activation gating variable of current Na
	hNa = x[3] # Inactivation gating variable of current Na
	mKd = x[4] # Activation gating variable of current Kd
	mCaL = x[5] # Activation gating variable of current CaL
	mCaN = x[6] # Activation gating variable of current CaN
	oERG = x[7] # ERG potassium current activation
	iERG = x[8] # ERG current intermediate
	gCaLi = x[9] # Intracellular maximum conductance of current CaL
	gCaL = x[10] # Maximum conductance of current CaL
	zCaL = x[11] # Integral variable of current CaL
	gCaNi = x[12] # Intracellular maximum conductance of current CaN
	gCaN = x[13] # Maximum conductance of current CaN
	zCaN = x[14] # Integral variable of current CaN

	# ODEs
	dx[1] = (1/C) * (- gNa*mNa^3*hNa^1*(V - 60.0) -
					   gKd*mKd^3*(V - -85.0) -
					   gCaL*mCaL^2*(V - 60.0) -
					   gCaN*mCaN^1*(V - 60.0) -
					   gERG*oERG^1*(V - -85.0) -
					   gNMDA*DA_NMDA_inf(V, Mg)^1*(V - 0.0) -
					   gleak*(V - -50.0) + Iapp)
	dx[2] = (1/DA_tau_mNa(V)) * (DA_mNa_inf(V) - mNa)
	dx[3] = (1/DA_tau_hNa(V)) * (DA_hNa_inf(V) - hNa)
	dx[4] = (1/DA_tau_mKd(V)) * (DA_mKd_inf(V) - mKd)
	dx[5] = (1/DA_tau_mCaL(V)) * (DA_mCaL_inf(V) - mCaL)
	dx[6] = (1/DA_tau_mCaN(V)) * (DA_mCaN_inf(V) - mCaN)
	dx[7] = a0ERG(V) * (1-oERG-iERG) + biERG(V) * iERG - oERG * (aiERG(V)+b0ERG(V))
    dx[8] = aiERG(V) * oERG - biERG(V) * iERG

	# Retrieving which line of the sensitivity matrix matter
	timescales = [2, 3]

	# Retrieving which column of the sensitivity matrix belong to unmodulated conductances
	unmodulated = [1, 2, 5, 6]

	# Retrieving which column of the sensitivity matrix belong to modulated conductances
	modulated = [3, 4]

	# Computing the right hand side of the linear system
	gDICr = [gsth, guth]
	gDICr = gDICr - SVth[timescales, unmodulated] * collect(p[2:5])

	# Computing the left hand side of the linear system
	Smod = SVth[timescales, modulated]

	# Computing the solution of the linear system
	g_r = \(Smod, gDICr)

	# Error signals and control inputs
	eCaL = g_r[1] - gCaL
	uCaL = Kp * eCaL + Ki * zCaL
	eCaN = g_r[2] - gCaN
	uCaN = Kp * eCaN + Ki * zCaN

	# ODEs of the controller
	dx[9] = α * gCaL - α * gCaLi - β * gCaLi + uCaL
	dx[10] = α * gCaLi - α * gCaL
	dx[11] = eCaL
	dx[12] = α * gCaN - α * gCaNi - β * gCaNi + uCaN
	dx[13] = α * gCaNi - α * gCaN
	dx[14] = eCaN
end
```

Note that a small adapation of `ControlledDA_ODE.jl` had been realized so that the name of the gating functions match the actual ones. Moreover, NMDA current had been changed to be instantaneaous (no gating variable dynamic) and the ERG gating variable dynamics were corrected since this current dynamic is non conventional.
