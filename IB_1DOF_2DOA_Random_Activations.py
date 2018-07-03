from pendulum_eqns.sim_eqns_ActIB_random_activations import *
from useful_functions import *

TotalX,TotalU = run_N_sim_rand_act(NumberOfTrials=1)
figs = plot_N_sim_rand_act(Time,TotalX,TotalU,Return=True)

plt.show()

# save_figures("1DOF_2DOA_Random_Activations")
# plt.close('all')
