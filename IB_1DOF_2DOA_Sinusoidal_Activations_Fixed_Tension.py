from pendulum_eqns.sim_eqns_ActIB_sinusoidal_activations import *
from useful_functions import *
import pickle
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
from scipy.stats import pearsonr
from danpy.useful_functions import is_number, save_figures

X_o = np.array([r(0),dr(0)])
InitialTensions = return_initial_tension(
                        X_o,
                        ReturnMultipleInitialTensions=True,
                        Bounds = [[0.15*F_MAX1,0.5*F_MAX1],[0.15*F_MAX2,0.5*F_MAX2]],
                        InitialAngularAcceleration=0
                        ) # len::10
InitialTensions=[InitialTensions[5]]
NumberOfTensionTrials = 1
NumberOfTrials = 25
InitialTensionsFromSuccessfulTrials = []
TerminalWidth = get_terminal_width()
count = 0
InitialMuscleLengths = []
for i in range(NumberOfTensionTrials):
    try:
        TensionTrialTitle = (
            "          Tension Setting "
            + str(i+1)
            + "/" +str(NumberOfTensionTrials)
            + "          \n")
        print(
        	" "*int(TerminalWidth/2 - len(TensionTrialTitle)/2)
        	+ colored(TensionTrialTitle,'blue',attrs=["underline","bold"])
        	)

        TotalX_temp,TotalU_temp = run_N_sim_IB_sinus_act(
                NumberOfTrials=NumberOfTrials,
                FixedInitialTension=InitialTensions[i],
                Amp="Scaled",
                Freq=1,
                PhaseOffset=-np.pi/3,
                InitialAngularAcceleration=0,
                InitialAngularSnap=0
                )
        _,Error_temp = plot_N_sim_IB_sinus_act(
                Time,TotalX_temp,
                TotalU_temp,Return=True,
                ReturnError=True)
        plt.close('all')
        count+=1

        if count == 1:
            TotalX = TotalX_temp
            TotalU = TotalU_temp
            Error1 = Error_temp[0]
            Error2 = Error_temp[1]
            Error = [Error1,Error2]
            InitialTensionsFromSuccessfulTrials.append(TotalX_temp[0,2:4,0])
            InitialMuscleLengths.append(TotalX_temp[0,4:6,0])
        else:
            TotalX = np.concatenate([TotalX,TotalX_temp],axis=0)
            TotalU = np.concatenate([TotalU,TotalU_temp],axis=0)
            Error1 = np.concatenate([Error1,Error_temp[0]],axis=0)
            Error2 = np.concatenate([Error2,Error_temp[1]],axis=0)
            Error = [Error1,Error2]
            InitialTensionsFromSuccessfulTrials.append(TotalX_temp[0,2:4,0])
            InitialMuscleLengths.append(TotalX_temp[0,4:6,0])
    except:
        print("Trial " + str(i+1) + " Failed...")

print("Number of Total Trials: " + str(NumberOfTensionTrials) + "\n")
print(
    "Number of Successful Trials: "
    + str(len(InitialTensionsFromSuccessfulTrials))
    )

if len(InitialTensions) != 0:
    figs = plot_N_sim_IB_sinus_act(Time,TotalX,TotalU,Return=True)

    additional_figs = plot_l_m_approximation_error_vs_tendon_tension(
                            Time,TotalX,
                            Error,Return=True,
                            InitialTensions=InitialTensionsFromSuccessfulTrials
                            )

    lm1o = np.array([TotalX[i,4,0] for i in range(TotalX.shape[0])])
    lm2o = np.array([TotalX[i,5,0] for i in range(TotalX.shape[0])])

    MAE1 = np.array([np.mean(abs(Error[0][i,:])) for i in range(Error[0].shape[0])])
    MAE2 = np.array([np.mean(abs(Error[1][i,:])) for i in range(Error[0].shape[0])])

    figMAE,(ax1,ax2) = plt.subplots(1,2,figsize=(10,5))
    plt.subplots_adjust(bottom=0.2)
    ax1.spines["top"].set_visible(False)
    ax1.spines['right'].set_visible(False)
    ax2.spines["top"].set_visible(False)
    ax2.spines['right'].set_visible(False)
    ax1.set_title("Muscle 1", fontsize=16)
    ax2.set_title("Muscle 2", fontsize=16)
    ax1.set_xlabel("Initial Normalized\nMuscle Fascicle Length",fontsize=14)
    ax1.set_ylabel("Percent Mean Absolute Error",fontsize=14)
    ax1.scatter(lm1o/lo1,100*(MAE1/lo1))
    ax1.text(
        0.5,0.9, f"PCC = {pearsonr(lm1o,MAE1)[0]:0.3f}",
        transform=ax1.transAxes,
        horizontalalignment='center',
        verticalalignment='center',
        color = "k",
        fontsize=14,
        bbox=dict(
            boxstyle='round,pad=0.5',
            edgecolor='k',
            facecolor='w'
        )
    )
    ax2.scatter(lm2o/lo2,100*(MAE2/lo2))
    ax2.text(
        0.5,0.9, f"PCC = {pearsonr(lm2o,MAE2)[0]:0.3f}",
        transform=ax2.transAxes,
        horizontalalignment='center',
        verticalalignment='center',
        color = "k",
        fontsize=14,
        bbox=dict(
            boxstyle='round,pad=0.5',
            edgecolor='k',
            facecolor='w'
        )
    )
    ax1.set_ylim([0,2])
    ax2.set_ylim([0,2])
    # plt.show()

    folderPath = save_figures(
        "output_figures/integrator_backstepping_sinusoidal_activations_fixed_tensions/",
        "1DOF_2DOA",
        {
            "Initial Muscle Lengths" : InitialMuscleLengths,
            "Initial Tendon Tensions" : InitialTensions
        },
        returnPath=True,
        saveAsPDF=True,
        saveAsMD=True
    )
    # save_figures("output_figures/integrator_backstepping_sinusoidal_activations_fixed_muscle_lengths/","1DOF_2DOA_v1.0",SaveAsPDF=True)
    plt.close('all')
    FormatedSaveData = {
            "States" : TotalX,
            "Input" : TotalU,
            "Error" : Error,
            "Initial Tensions" : InitialTensionsFromSuccessfulTrials
            }
    pickle.dump(
        FormatedSaveData,
        open(
            folderPath/"output.pkl",
            "wb"
            )
        )
else:
    print("All Trials Unsuccessful...")
