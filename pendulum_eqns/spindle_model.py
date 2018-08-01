import numpy as np
from sim_eqns_ActIB_gaussian_activations_around_previous_input import *

def return_β(β0,β1,β2,f_dynamic,f_static):
    result =  β0 + β1*f_dynamic + β2*f_static
    return(result)

def return_C(NormalizedFiberVelocity):
    C_L = 1
    C_S = 0.4200
    result = (
        (NormalizedFiberVelocity>=0)*C_L
        + (NormalizedFiberVelocity<0)*C_S
        )
    return(result)

GammaGainMatrix = np.round((10+np.random.rand(2,2)*(250-10))/5)*5

γ_dynamic_1 = GammaGainMatrix[0,0]
γ_static_1 = GammaGainMatrix[1,0]

γ_dynamic_2 = GammaGainMatrix[0,1]
γ_static_2 = GammaGainMatrix[1,1]

K_sr = 10.4649
K_pr = 0.1500
M = 0.0002 # in (Force Units)/(L_o/s²)
L_sr_threshold = 0.0423
L_pr_threshold = 0.89
a = 0.3
R = 0.46
L_sr_o = 0.04
L_pr_o = 0.76
p = 2

def update_f_dynamic_bag_1(i):
    """
    γ_dynamic_1 is a fixed value

    f_dynamic_bag_1 is an array

    i is the current timestep

    t is the time array (for generating dt)
    """

    dt = t[1]-t[0]
    tau = 0.149
    Freq_bag1 = 60
    df_dynamic_bag_1 = (γ_dynamic_1**p/(γ_dynamic_1**p + Freq_bag1**2) - f_dynamic_bag_1[i])/tau
    f_dynamic_bag_1[i+1] = f_dynamic_bag_1[i] + df_dynamic_bag_1*dt

def return_f_dynamic_bag_2(
        i,
        t,
        γ_static,
        f_static
        ):
    """
    γ_static is a fixed value

    f_static is an array

    i is the current timestep

    t is the time array (for generating dt)
    """

    dt = t[1]-t[0]
    tau = 0.149
    Freq_bag2 = 60
    df_static = (γ_static**p/(γ_static**p + Freq_bag2**2) - f_static[i])/tau
    f_static[i+1] = f_static[i] + df_static*dt
    return(f_static)

def return_f_dynamic_chain(
        i,
        γ_static,
        f_static
        ):
    """
    γ_static is a fixed value

    f_static is an array

    i is the current timestep
    """

    Freq_chain = 90
    f_static = (γ_static**p/(γ_static**p + Freq_chain**2)
    return(f_static)

def find_initial_fiber_tension_bag_1(
        f_dynamic_o,
        InitialNormalizedFiberLength):

    Gamma = Gamma1*f_dynamic_o
    InitialTension = (
        (K_pr*K_sr/(K_pr + K_sr))
        * (InitialNormalizedFiberLength - L_sr_o - L_pr_o + Gamma/K_pr)
        )
    return(InitialTension)

def update_intrafusal_fiber_tension_bag_1(
        i,
        NormalizedFiberLength,
        NormalizedFiberVelocity,
        NormalizedFiberAcceleration):

    C = return_C(NormalizedFiberVelocity[i])
    β = return_β(0.0605,0.2592,0,f_dynamic_bag_1[i],_)

    update_f_dynamic_bag_1(i)
    Gamma = Gamma1*f_dynamic_bag_1[i+1]

    FiberTensionAcceleration_bag_1[i+1] = \
            (K_sr/M)
                * ( C
                        * β
                        * np.sign(NormalizedFiberVelocity[i] - FiberTensionVelocity_bag_1[i]/K_sr)
                        * (abs(NormalizedFiberVelocity[i] - FiberTensionVelocity_bag_1[i]/K_sr)**a)
                        * (NormalizedFiberLength[i] - L_sr_o - FiberTension_bag_1[i]/K_sr - R)
                    + K_pr*(NormalizedFiberLength[i] - L_sr_o - FiberTension_bag_1[i]/K_sr - L_pr_o)
                    + M*NormalizedFiberAcceleration[i]
                    + Gamma
                    - FiberTension_bag_1[i]
                )
    FiberTensionVelocity_bag_1[i+1] = FiberTensionVelocity_bag_1[i] + FiberTensionAcceleration_bag_1[i+1]*dt
    FiberTension_bag_1[i+1] = FiberTension_bag_1[i] + FiberTensionVelocity_bag_1[i+1]*dt

def update_primary_afferent(i,FiberTension):
    """
    FiberTension should be FiberTension_bag_1[i+1]
    """
    Af_primary[i+1] = G*(FiberTension/K_sr - (L_sr_threshold-L_sr_o))

def run_this_shit(
        NormalizedFiberLength,
        NormalizedFiberVelocity,
        NormalizedFiberAcceleration
        ):

    f_dynamic_bag_1 = np.zeros(1,len(Time))
    f_static_bag_2 = np.zeros(1,len(Time))
    f_static_chain = np.zeros(1,len(Time))
    FiberTension_bag_1 = np.zeros(1,len(Time))
    FiberTension_bag_1[0] = find_initial_fiber_tension_bag_1(
                                        NormalizedFiberLength[0],
                                        f_dynamic_bag_1[0]
                                        )
    FiberTensionVelocity_bag_1 = np.zeros(1,len(Time))
    FiberTensionAcceleration_bag_1 = np.zeros(1,len(Time))
    Af_primary = np.zeros(1,len(Time))

    for i in range(len(Time)):
        update_intrafusal_fiber_tension_bag_1(
                i,
                NormalizedFiberLength,
                NormalizedFiberVelocity,
                NormalizedFiberAcceleration)
        update_primary_afferent(i,FiberTension_bag_1[i+1])

    return(Af_primary)
