import numpy as np
import matplotlib.pyplot as plt
from danpy.sb import dsb
# from scipy import integrate
from pendulum_eqns.state_equations import *

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

# GammaGainMatrix = np.round((10+np.random.rand(2,2)*(250-10))/5)*5
#
# gamma_dynamic_1 = GammaGainMatrix[0,0]
# gamma_static_1 = GammaGainMatrix[1,0]
#
# gamma_dynamic_2 = GammaGainMatrix[0,1]
# gamma_static_2 = GammaGainMatrix[1,1]

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

Gamma1 = 0.0289
#
# def return_f_dynamic_bag_2(
#         i,
#         t,
#         gamma_static,
#         f_static
#         ):
#     """
#     gamma_static is a fixed value
#
#     f_static is an array
#
#     i is the current timestep
#
#     t is the time array (for generating dt)
#     """
#
#     dt = t[1]-t[0]
#     tau = 0.149
#     Freq_bag2 = 60
#     df_static = (gamma_static**p/(gamma_static**p + Freq_bag2**2) - f_static[i])/tau
#     f_static[i+1] = f_static[i] + df_static*dt
#     return(f_static)
#
# def return_f_dynamic_chain(
#         i,
#         gamma_static,
#         f_static
#         ):
#     """
#     gamma_static is a fixed value
#
#     f_static is an array
#
#     i is the current timestep
#     """
#
#     Freq_chain = 90
#     f_static = (gamma_static**p/(gamma_static**p + Freq_chain**2))
#     return(f_static)

def find_initial_fiber_tension_bag_1(
        f_dynamic_o,
        InitialNormalizedFiberLength):

    Gamma = Gamma1*f_dynamic_o
    InitialTension = (
        (K_pr*K_sr/(K_pr + K_sr))
        * (InitialNormalizedFiberLength - L_sr_o - L_pr_o + Gamma/K_pr)
        )
    return(InitialTension)

gamma_dynamic_1 = np.round((10+np.random.rand()*(250-10))/5)*5

def return_primary_afferent_from_single_trial(
        Time,
        NormalizedFiberLength,
        NormalizedFiberVelocity,
        gamma_dynamic = 100
        ):

    gamma_dynamic_1 = gamma_dynamic
    def return_f_dynamic_bag_1_function(gamma_dynamic):
        tau = 0.149
        Freq_bag1 = 60
        def f_dynamic(t):
            return((1-np.exp(-t/0.149))*(gamma_dynamic**p/(gamma_dynamic**p + Freq_bag1**p)))
        return(f_dynamic)
    f_dynamic_bag_1 = return_f_dynamic_bag_1_function(gamma_dynamic)(Time)
    # def update_f_dynamic_bag_1(i):
    #     """
    #     gamma_dynamic_1 is a fixed value
    #
    #     f_dynamic_bag_1 is an array
    #
    #     i is the current timestep
    #
    #     t is the time array (for generating dt)
    #     """
    #
    #     dt = Time[1]-Time[0]
    #     tau = 0.149
    #     Freq_bag1 = 60
    #     df_dynamic_bag_1 = (gamma_dynamic_1**p/(gamma_dynamic_1**p + Freq_bag1**p) - f_dynamic_bag_1[i])/tau
    #     f_dynamic_bag_1[i+1] = f_dynamic_bag_1[i] + df_dynamic_bag_1*dt
    #

    def update_intrafusal_fiber_tension_bag_1(
            i,
            NormalizedFiberLength,
            NormalizedFiberVelocity,
            NormalizedFiberAcceleration):

        C = return_C(NormalizedFiberVelocity[i])
        β = return_β(0.0605,0.2592,0,f_dynamic_bag_1[i],0)

        # update_f_dynamic_bag_1(i)
        Gamma = Gamma1*f_dynamic_bag_1[i]
        # test =  (C
        #         * β
        #         * np.sign(NormalizedFiberVelocity[i] - FiberTensionVelocity_bag_1[i]/K_sr)
        #         * (abs(NormalizedFiberVelocity[i] - FiberTensionVelocity_bag_1[i]/K_sr)**a))
        FiberTensionAcceleration_bag_1[i+1] = \
                (K_sr/M) \
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
        # return(test)

    def update_primary_afferent(i,FiberTension):
        """
        FiberTension should be FiberTension_bag_1[i+1]
        """
        G = 20000
        Af_primary[i+1] = G*(FiberTension/K_sr - (L_sr_threshold-L_sr_o))

    dt = Time[1]-Time[0]
    NormalizedFiberAcceleration = np.gradient(NormalizedFiberVelocity,dt)

    # f_dynamic_bag_1 = np.zeros(len(Time))
    # f_static_bag_2 = np.zeros(len(Time))
    # f_static_chain = np.zeros(len(Time))
    FiberTension_bag_1 = np.zeros(len(Time))
    FiberTension_bag_1[0] = find_initial_fiber_tension_bag_1(
                                        NormalizedFiberLength[0],
                                        f_dynamic_bag_1[0]
                                        )
    FiberTensionVelocity_bag_1 = np.zeros(len(Time))
    FiberTensionAcceleration_bag_1 = np.zeros(len(Time))
    Af_primary = np.zeros(len(Time))
    # TestVar = np.zeros(len(Time))

    for i in range(len(Time)-1):
        update_intrafusal_fiber_tension_bag_1(
                i,
                NormalizedFiberLength,
                NormalizedFiberVelocity,
                NormalizedFiberAcceleration)
        # TestVar[i]=test
        update_primary_afferent(i,FiberTension_bag_1[i+1])

    return(Af_primary)


def return_Ia_error(
        N,
        Time,
        TotalX,
        **kwargs
        ):

    NormalizeError = kwargs.get("NormalizeError",False)
    assert type(NormalizeError) == bool, "NormalizeError must be either True or False"

    ReturnAllData = kwargs.get("ReturnAllData",False)
    assert type(ReturnAllData) == bool, "ReturnAllData must be either True or False"

    fig1 = plt.figure()
    ax1 = plt.gca()
    ax1.set_title("Muscle 1")
    ax1.set_xlabel(r"$\gamma_{dynamic}$")
    if NormalizeError == True:
        ax1.set_ylabel("Normalized Average Error Value (pps)")
    else:
        ax1.set_ylabel("Average Error Value (pps)")

    fig2 = plt.figure()
    ax2 = plt.gca()
    ax2.set_title("Muscle 2")
    ax2.set_xlabel(r"$\gamma_{dynamic}$")
    if NormalizeError == True:
        ax2.set_ylabel("Normalized Average Error Value (pps)")
    else:
        ax2.set_ylabel("Average Error Value (pps)")

    # ApproxTestVarArray1 = np.zeros((np.shape(TotalX)[0],N,len(Time)))
    # TestVarArray1 = np.zeros((np.shape(TotalX)[0],N,len(Time)))

    # ApproxTestVarArray2 = np.zeros((np.shape(TotalX)[0],N,len(Time)))
    # TestVarArray2 = np.zeros((np.shape(TotalX)[0],N,len(Time)))

    AllPrimaryAfferent1 = np.zeros((np.shape(TotalX)[0],N,len(Time)))
    AllPrimaryAfferent2 = np.zeros((np.shape(TotalX)[0],N,len(Time)))
    AllApproxPrimaryAfferent1 = np.zeros((np.shape(TotalX)[0],N,len(Time)))
    AllApproxPrimaryAfferent2 = np.zeros((np.shape(TotalX)[0],N,len(Time)))

    for j in range(np.shape(TotalX)[0]):
        L_approx_1 = (integrate.cumtrapz(
                         np.array(list(map(lambda X: v_MTU1(X),TotalX[j,:,:].T))),
                         Time,
                         initial=0
                         )
                    + np.ones(len(Time))*TotalX[j,4,0]
                    )/lo1
        dL_approx_1 = np.array(list(map(lambda X: v_MTU1(X),TotalX[j,:,:].T)))/lo1

        GammaArray_1 = np.zeros(N)
        AverageErrorArray_1 = np.zeros(N)
        gamma_dynamic_array = np.linspace(10,250,N)
        statusbar1 = dsb(0,N,title="Muscle 1, Trial " + str(j+1) + "/" + str(np.shape(TotalX)[0]))
        for i in range(N):
            gamma_dynamic_1 = gamma_dynamic_array[i]
            GammaArray_1[i] = gamma_dynamic_1
            Approx_Aff_potential_bag_1 = return_primary_afferent_from_single_trial(Time,L_approx_1,dL_approx_1,gamma_dynamic=gamma_dynamic_1)
            Aff_potential_bag_1 = return_primary_afferent_from_single_trial(Time,TotalX[j,4,:]/lo1,TotalX[j,6,:]/lo1,gamma_dynamic=gamma_dynamic_1)

            AllPrimaryAfferent1[j,i,:] = Aff_potential_bag_1
            AllApproxPrimaryAfferent1[j,i,:] = Approx_Aff_potential_bag_1

            AverageError1 = \
                np.average(
                    abs(
                        Aff_potential_bag_1[int(len(Time)*0.1):]
                        - Approx_Aff_potential_bag_1[int(len(Time)*0.1):]
                    )
                )
            if NormalizeError == True:
                PeakToPeak = np.array(
                    [max(Aff_potential_bag_1[i:i+10000])-min(Aff_potential_bag_1[i:i+10000]) for i in range(100,len(Time)-10000,100)]
                    ).mean()
                AverageError1 = AverageError1/PeakToPeak
            AverageErrorArray_1[i] = AverageError1
            statusbar1.update(i)

        L_approx_2 = (integrate.cumtrapz(
                         np.array(list(map(lambda X: v_MTU2(X),TotalX[j,:,:].T))),
                         Time,
                         initial=0
                         )
                    + np.ones(len(Time))*TotalX[j,5,0]
                    )/lo2
        dL_approx_2 = np.array(list(map(lambda X: v_MTU2(X),TotalX[j,:,:].T)))/lo2

        GammaArray_2 = np.zeros(N)
        AverageErrorArray_2 = np.zeros(N)
        statusbar2 = dsb(0,N,title="Muscle 2, Trial " + str(j+1) + "/" + str(np.shape(TotalX)[0]))
        for i in range(N):
            gamma_dynamic_1 = gamma_dynamic_array[i]
            # gamma_dynamic_1 = np.round((10+np.random.rand()*(250-10))/5)*5
            GammaArray_2[i] = gamma_dynamic_1
            Approx_Aff_potential_bag_1 = return_primary_afferent_from_single_trial(Time,L_approx_2,dL_approx_2,gamma_dynamic=gamma_dynamic_1)
            Aff_potential_bag_1 =  return_primary_afferent_from_single_trial(Time,TotalX[j,5,:]/lo2,TotalX[j,7,:]/lo2,gamma_dynamic=gamma_dynamic_1)

            AllPrimaryAfferent2[j,i,:] = Aff_potential_bag_1
            AllApproxPrimaryAfferent2[j,i,:] = Approx_Aff_potential_bag_1

            AverageError2 = \
                np.average(
                    abs(
                        Aff_potential_bag_1[int(len(Time)*0.1):]
                        - Approx_Aff_potential_bag_1[int(len(Time)*0.1):]
                    )
                )
            if NormalizeError == True:
                PeakToPeak = np.array(
                    [max(Aff_potential_bag_1[i:i+10000])-min(Aff_potential_bag_1[i:i+10000]) for i in range(100,len(Time)-10000,100)]
                    ).mean()
                AverageError2 = AverageError2/PeakToPeak
            AverageErrorArray_2[i] = AverageError2
            statusbar2.update(i)

        # plt.figure()
        ax1.plot(GammaArray_1,AverageErrorArray_1)
        # ax1 = plt.gca()


        # plt.figure()
        ax2.plot(GammaArray_2,AverageErrorArray_2)
        # ax2 = plt.gca()

    plt.show()

    if ReturnAllData == True:
        Data = {
            "Error is Normalized?" : NormalizeError,
            "Muscle 1" : {
                "Gamma Dynamic Gains" : GammaArray_1,
                "Primary Afferent Signals" : AllPrimaryAfferent1,
                "Approximate Primary Afferent Signals" : AllApproxPrimaryAfferent1
            },
            "Muscle 2" : {
                "Gamma Dynamic Gains" : GammaArray_2,
                "Primary Afferent Signals" : AllPrimaryAfferent2,
                "Approximate Primary Afferent Signals" : AllApproxPrimaryAfferent2
            },
        }
        return(Data)
