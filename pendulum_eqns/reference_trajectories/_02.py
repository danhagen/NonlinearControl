import numpy as np
"""
Continuously differentiable up to d2r
"""
Amp = 7.5*np.pi/180
Base = 90*np.pi/180
Freq = 2*np.pi
Delay1 = 5*np.pi/(2*Freq)
# N_seconds = 1
# N = N_seconds*10000 + 1
# t = np.linspace(0,N_seconds,N)
# dt = t[1]-t[0]

### Reference Trajectory ###

def r(t):
    if 0<=t<Delay1:
        return(np.pi*Amp/4 + Base)
    elif Delay1<=t<(Delay1+(1/4)*(2*np.pi/Freq)):
        return((Amp/4)*np.sin(2*Freq*(t-Delay1)) - Amp*Freq*t/2 + 3*np.pi*Amp/2 + Base)
    else:
        return(-Amp*np.sin(Freq*(t-Delay1-(1/4)*(2*np.pi/Freq))) + Base)
def dr(t):
    if 0<=t<Delay1:
        return(0)
    elif Delay1<=t<(Delay1+(1/4)*(2*np.pi/Freq)):
        return((Amp*Freq/2)*np.cos(2*Freq*(t-Delay1)) - Amp*Freq/2)
    else:
        return(-Amp*Freq*np.cos(Freq*(t-Delay1-(1/4)*(2*np.pi/Freq))))
def d2r(t):
    if 0<=t<Delay1:
        return(0)
    elif Delay1<=t<(Delay1+(1/4)*(2*np.pi/Freq)):
        return(-(Amp*Freq**2)*np.sin(2*Freq*(t-Delay1)))
    else:
        return(Amp*Freq**2*np.sin(Freq*(t-Delay1-(1/4)*(2*np.pi/Freq))))
def d3r(t):
    if 0<=t<Delay1:
        return(0)
    elif Delay1<=t<(Delay1+(1/4)*(2*np.pi/Freq)):
        return(-(2*Amp*Freq**3)*np.cos(2*Freq*(t-Delay1)))
    else:
        return(Amp*Freq**3*np.cos(Freq*(t-Delay1-(1/4)*(2*np.pi/Freq))))
def d4r(t):
    if 0<=t<Delay1:
        return(0)
    elif Delay1<=t<(Delay1+(1/4)*(2*np.pi/Freq)):
        return((4*Amp*Freq**4)*np.sin(2*Freq*(t-Delay1)))
    else:
        return(-Amp*Freq**4*np.sin(Freq*(t-Delay1-(1/4)*(2*np.pi/Freq))))

############################
