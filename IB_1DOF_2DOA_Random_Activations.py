from pendulum_eqns.sim_eqns_ActIB_random_activations import *
from useful_functions import *

TotalX,TotalU = run_N_sim_rand_act(NumberOfTrials=3)
figs = plot_N_sim_rand_act(Time,TotalX,TotalU,Return=True)

# plt.show()

save_figures("output_figures/random_activations/","1DOF_2DOA_v1.0",SaveAsPDF=True)
plt.close('all')
