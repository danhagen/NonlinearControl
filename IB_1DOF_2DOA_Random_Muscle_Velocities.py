from pendulum_eqns.sim_eqns_VmIB_random_muscle_velocities import *
from useful_functions import *

TotalX,TotalU = run_N_sim_rand_Vm(NumberOfTrials=1)
figs = plot_N_sim_rand_Vm(Time,TotalX,TotalU,Return=True)

plt.show()

# save_figures("1DOF_2DOA_Random_Activations")
# plt.close('all')
