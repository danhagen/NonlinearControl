from pendulum_eqns.sim_eqns_ActIB_gaussian_activations_around_previous_input import *
from useful_functions import *

X_o = np.array([r(0),dr(0)])
InitialTension = return_initial_tension(X_o)

TotalX,TotalU = run_N_sim_gauss_act(NumberOfTrials=100,FixedInitialTension=InitialTension)
figs = plot_N_sim_gauss_act(Time,TotalX,TotalU,Return=True)

# plt.show()

save_figures("output_figures/fixed_initial_tension_gauss_act/" + "1DOF_2DOA_v1.0")
plt.close('all')
