from pendulum_eqns.sim_eqns_ActIB_gaussian_activations_around_previous_input import *
from useful_functions import *

TotalX,TotalU = run_N_sim_gauss_act(NumberOfTrials=3)
figs = plot_N_sim_gauss_act(Time,TotalX,TotalU,Return=True)

# plt.show()

save_figures("output_figures/gaussian_activations/","1DOF_2DOA_v1.0",SaveAsPDF=True)
plt.close('all')
