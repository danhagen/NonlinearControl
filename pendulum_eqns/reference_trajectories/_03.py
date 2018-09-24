import numpy as np
"""
Continuously differentiable up to d2r
"""
Amp = 7.5*np.pi/180
Base = 90*np.pi/180
InitialAngle = Base-Amp
FinalAngle = Base+Amp
Freq = 2*np.pi
Period = 2*np.pi/Freq
HalfPeriod = Period/2
d_tau_dt = 1/HalfPeriod
# N_seconds = 1
# N = N_seconds*10000 + 1
# t = np.linspace(0,N_seconds,N)
# dt = t[1]-t[0]

### Reference Trajectory ###

def r(t):
    tau = (t - np.floor(t*d_tau_dt)*HalfPeriod)/HalfPeriod
    if int(np.floor(t*d_tau_dt))%2 == 0:
        return(InitialAngle + (FinalAngle-InitialAngle)*(6*tau**5 - 15*tau**4 + 10*tau**3))
    else:
        return(FinalAngle + (InitialAngle-FinalAngle)*(6*tau**5 - 15*tau**4 + 10*tau**3))
def dr(t):
    tau = (t - np.floor(t*d_tau_dt)*HalfPeriod)/HalfPeriod
    if int(np.floor(t*d_tau_dt))%2 == 0:
        return((FinalAngle-InitialAngle)*(30*tau**4 -60*tau**3 + 30*tau**2))
    else:
        return((InitialAngle-FinalAngle)*(30*tau**4 -60*tau**3 + 30*tau**2))
def d2r(t):
    tau = (t - np.floor(t*d_tau_dt)*HalfPeriod)/HalfPeriod
    if int(np.floor(t*d_tau_dt))%2 == 0:
        return((FinalAngle-InitialAngle)*(120*tau**3 - 180*tau**2 + 60*tau))
    else:
        return((InitialAngle-FinalAngle)*(120*tau**3 - 180*tau**2 + 60*tau))
def d3r(t):
    tau = (t - np.floor(t*d_tau_dt)*HalfPeriod)/HalfPeriod
    if int(np.floor(t*d_tau_dt))%2 == 0:
        return((FinalAngle-InitialAngle)*(360*tau**2 - 360*tau + 60))
    else:
        return((InitialAngle-FinalAngle)*(360*tau**2 - 360*tau + 60))
def d4r(t):
    tau = (t - np.floor(t*d_tau_dt)*HalfPeriod)/HalfPeriod
    if int(np.floor(t*d_tau_dt))%2 == 0:
        return((FinalAngle-InitialAngle)*(720*tau - 360))
    else:
        return((InitialAngle-FinalAngle)*(720*tau - 360))

############################
