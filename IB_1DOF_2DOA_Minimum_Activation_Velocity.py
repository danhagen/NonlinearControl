from pendulum_eqns.sim_eqns_ActIB_minimize_activation_velocity import *
from useful_functions import *

TotalX,TotalU = run_N_sim_MAV(NumberOfTrials=1)
figs = plot_N_sim_MAV(Time,TotalX,TotalU,Return=True)

plt.show()

# save_figures("1DOF_2DOA_Minimum_Activation_Velocity")
# plt.close('all')
