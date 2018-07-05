from pendulum_eqns.sim_eqns_ActIB_minimize_activation_velocity import *
from useful_functions import *

TotalX,TotalU = run_N_sim_MAV(NumberOfTrials=3)
figs = plot_N_sim_MAV(Time,TotalX,TotalU,Return=True)

# plt.show()

save_figures("output_figures/minimum_activation_velocity/" + "1DOF_2DOA_v1.0")
plt.close('all')
