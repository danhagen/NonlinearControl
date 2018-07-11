from pendulum_eqns.sim_eqns_ActIB_gaussian_activations_around_previous_input import *
from useful_functions import *

NumberOfTensionTrials = 4
InitialTensions = []
for i in range(NumberOfTensionTrials):
    try:
        X_o = np.array([r(0),dr(0)])
        InitialTension = return_initial_tension(X_o)

        TotalX_temp,TotalU_temp = run_N_sim_gauss_act(NumberOfTrials=3,FixedInitialTension=InitialTension)
        _,Error_temp = plot_N_sim_gauss_act(
                Time,TotalX_temp,
                TotalU_temp,Return=True,
                ReturnError=True)
        plt.close('all')

        if i == 0:
            TotalX = TotalX_temp
            TotalU = TotalU_temp
            Error1 = Error_temp[0]
            Error2 = Error_temp[1]
            Error = [Error1,Error2]
            InitialTensions.append(TotalX_temp[0,2:4,0])
        else:
            TotalX = np.concatenate([TotalX,TotalX_temp],axis=0)
            TotalU = np.concatenate([TotalU,TotalU_temp],axis=0)
            Error1 = np.concatenate([Error1,Error_temp[0]],axis=0)
            Error2 = np.concatenate([Error2,Error_temp[1]],axis=0)
            Error = [Error1,Error2]
            InitialTensions.append(TotalX_temp[0,2:4,0])
    except:
        print("Trial " + str(i+1) + " Failed...")

print("Number of Total Trials: " + str(NumberOfTensionTrials) + "\n")
print("Number of Successful Trials: " + str(len(InitialTensions)))

if len(InitialTensions) != 0:
    figs = plot_N_sim_gauss_act(Time,TotalX,TotalU,Return=True)

    additional_figs = plot_l_m_approximation_error_vs_tendon_tension(
                            Time,TotalX,
                            Error,Return=True,
                            InitialTensions=InitialTensions
                            )
    # plt.show()

    save_figures("output_figures/gauss_act_initial_tension_sweep/","1DOF_2DOA_v1.0",SaveAsPDF=True)
    plt.close('all')
else:
    print("All Trials Unsuccessful...")
