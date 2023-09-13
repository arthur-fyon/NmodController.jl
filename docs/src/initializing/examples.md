# Example of existing models in the litterature

Through this documentation, two conductance based model coming from the litterature wil be used with *NmodController.jl*: a stomatogastric (STG) neuron model (Liu et al., 1998 "A model neuron with activity-dependent conductances regulated by multiple calcium sensors") and an adapted dopaminergic (DA) neuron model (Drion et al., 2011 "How modeling can reconcile apparently discrepant experimental results: the case of pacemaking in dopaminergic neurons").

## Initializing the STG model
The STG model is composed of 7 voltage gated ionic currents which one is calcium dependent:
1. Transient sodium current $I_\mathrm{Na}$ (2 gating variables);
2. T-type calcium current $I_\mathrm{CaT}$ (2 gating variables);
3. Slow calcium current $I_\mathrm{CaS}$ (2 gating variables);
4. A-type potassium current $I_\mathrm{A}$ (2 gating variables);
5. Calcium controlled potassium current $I_\mathrm{KCa}$ (1 gating variable caclium dpendent);
6. Delayed rectified potassium current $I_\mathrm{Kd}$ (1 gating variable);
7. H current $I_\mathrm{H}$ (1 gating variable).

The voltage equation writes

$C \dot V = - \bar{g}_\mathrm{Na}m^3_\mathrm{Na}h_\mathrm{Na}(V-E_\mathrm{Na}) - \bar{g}_\mathrm{CaT}m^3_\mathrm{CaT}h_\mathrm{CaT}(V-E_\mathrm{Ca}) - \bar{g}_\mathrm{CaS}m^3_\mathrm{CaS}h_\mathrm{CaS}(V-E_\mathrm{Ca}) - \bar{g}_\mathrm{A}m^3_\mathrm{A}h_\mathrm{A}(V-E_\mathrm{K}) - \bar{g}_\mathrm{KCa}m^4_\mathrm{KCa}(V-E_\mathrm{K}) - \bar{g}_\mathrm{Kd}m^4_\mathrm{Kd}(V-E_\mathrm{K}) - \bar{g}_\mathrm{H}m_\mathrm{H}(V-E_\mathrm{H}) - g_\mathrm{leak}(V-E_\mathrm{leak}) + I_{ext}(t)$

and the intracellular calcium dynamic writes

$\tau_{Ca} \cdot \dot{\[Ca\]} = -0.94\cdot  I_{CaT} -0.94\cdot  I_{CaS} - \[Ca\] + Ca_\mathrm{eq}$

with $\tau_{Ca} = 20$ and $Ca_\mathrm{eq} = 0.05$.

The next few lines of code show how to initialize such model using *NmodController.jl*.