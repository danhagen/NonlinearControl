from IB_minimize_activation_velocity import *

TotalX,TotalU = run_N_sim_MAT(NumberOfTrials=1)
figs = plot_N_sim_MAT(t,TotalX,TotalU,Return=True)

plt.show()

# save_figures("1DOF_2DOA_Minimum_Activation_Velocity")
# plt.close('all')
