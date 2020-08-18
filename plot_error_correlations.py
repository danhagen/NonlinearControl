import pickle
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import pearsonr
from pendulum_eqns.physiology.muscle_params_BIC_TRI import *

#### Fixed Initial Tendon Tension

out=pickle.load( open( Path(r"C:\Users\hagen\Documents\Github\NonlinearControl\output_figures\integrator_backstepping_sinusoidal_activations_fixed_tensions\2020_05_23_112705\output.pkl"), "rb" ) )

Error = out["Error"]
States = out["States"]

lm1o = np.array([States[i,4,0] for i in range(States.shape[0])])
lm2o = np.array([States[i,5,0] for i in range(States.shape[0])])

MAE1 = np.array([np.mean(abs(Error[0][i,:])) for i in range(Error[0].shape[0])])
MAE2 = np.array([np.mean(abs(Error[1][i,:])) for i in range(Error[0].shape[0])])

fig1,(ax1,ax2) = plt.subplots(1,2,figsize=(10,5))
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

### Fixed Initial Muscle Length

out=pickle.load( open( Path(r"C:\Users\hagen\Documents\Github\NonlinearControl\output_figures\integrator_backstepping_sinusoidal_activations_fixed_muscle_lengths\2020_05_23_115050\output.pkl"), "rb" ) )

Error = out["Error"]
States = out["States"]

fT1o = np.array([States[i,2,0] for i in range(States.shape[0])])
fT2o = np.array([States[i,3,0] for i in range(States.shape[0])])

MAE1 = np.array([np.mean(abs(Error[0][i,:])) for i in range(Error[0].shape[0])])
MAE2 = np.array([np.mean(abs(Error[1][i,:])) for i in range(Error[0].shape[0])])

fig2,(ax1,ax2) = plt.subplots(1,2,figsize=(10,5))
plt.subplots_adjust(bottom=0.2)
ax1.spines["top"].set_visible(False)
ax1.spines['right'].set_visible(False)
ax2.spines["top"].set_visible(False)
ax2.spines['right'].set_visible(False)
ax1.set_title("Muscle 1", fontsize=16)
ax2.set_title("Muscle 2", fontsize=16)
ax1.set_xlabel("Initial Normalized\nTendon Force",fontsize=14)
ax1.set_ylabel("Percent Mean Absolute Error",fontsize=14)
ax1.scatter(fT1o/F_MAX1,100*(MAE1/lo1))
ax1.text(
    0.5,0.9, f"PCC = {pearsonr(fT1o,MAE1)[0]:0.3f}",
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
ax2.scatter(fT2o/F_MAX2,100*(MAE2/lo2))
ax2.text(
    0.5,0.9, f"PCC = {pearsonr(fT2o,MAE2)[0]:0.3f}",
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

plt.show()
