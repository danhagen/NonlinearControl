from pendulum_eqns.sim_eqns_TTIB_random_tensions import *
from useful_functions import *

TotalX,TotalU = run_N_sim_rand_TT(NumberOfTrials=3)
figs = plot_N_sim_rand_TT(Time,TotalX,TotalU,Return=True)

# plt.show()

save_figures("output_figures/random_tendon_tensions/" + "1DOF_2DOA_v1.0")
plt.close('all')
