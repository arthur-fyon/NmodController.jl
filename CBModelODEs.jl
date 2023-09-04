#=
	This file contains differential equations describing the STG model
=#
# Function that outputs values of variables derivatives
function CB_ODE(dx, x, p, t)
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
	V = x[1](t)
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
	dx[2] = (1/tau_mNa(V)) * (mNa_inf(V) - mNa)
	dx[3] = (1/tau_hNa(V)) * (hNa_inf(V) - hNa)
	dx[4] = (1/tau_mCaT(V)) * (mCaT_inf(V) - mCaT)
	dx[5] = (1/tau_hCaT(V)) * (hCaT_inf(V) - hCaT)
	dx[6] = (1/tau_mCaS(V)) * (mCaS_inf(V) - mCaS)
	dx[7] = (1/tau_hCaS(V)) * (hCaS_inf(V) - hCaS)
	dx[8] = (1/tau_mA(V)) * (mA_inf(V) - mA)
	dx[9] = (1/tau_hA(V)) * (hA_inf(V) - hA)
	dx[10] = (1/tau_mKCa(V)) * (mKCa_inf(V, Ca) - mKCa)
	dx[11] = (1/tau_mKd(V)) * (mKd_inf(V) - mKd)
	dx[12] = (1/tau_mH(V)) * (mH_inf(V) - mH)
	dx[13] = (-0.94*gCaT*mCaT^3*hCaT^1*(V - 80.0) + -0.94*gCaS*mCaS^3*hCaS^1*(V - 80.0) + -Ca + 0.05) / 20.0
end
