from pendulum_eqns.sim_eqns_VmIB_random_muscle_velocities import *
from useful_functions import *

TotalX,TotalU = run_N_sim_rand_Vm(NumberOfTrials=3)
figs = plot_N_sim_rand_Vm(Time,TotalX,TotalU,Return=True)

# plt.show()

save_figures("output_figures/random_muscle_velocities/","1DOF_2DOA_v1.0",SaveAsPDF=True)
plt.close('all')
