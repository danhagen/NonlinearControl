from pendulum_eqns.physiology.muscle_params_BIC_TRI import *
from pendulum_eqns.state_equations import *
from scipy.integrate import cumtrapz
import matplotlib.pyplot as plt
from danpy.sb import dsb

Theta_i = np.pi/6
Theta_f = 2*np.pi/3
Omega = 1
T_end = (Theta_f-Theta_i)/Omega
N = int(T_end*1000+1)
Time = np.linspace(0,T_end,N)
X = np.zeros((2,len(Time)))
X[0,:] = Theta_i + Omega*Time
X[1,:] = Omega*np.ones(np.shape(Time))

lTo1 = 0.5

def error_func_1(X,T,T_i,Time):
    assert np.shape(X)[0]>=2, "X must be an array of shape (M,N) where M>=2."
    error = (lTo1*kT/np.cos(α1))*np.log((np.exp(T/(F_MAX1*cT*kT))-1)/(np.exp(T_i/(F_MAX1*cT*kT))-1)) + ((np.cos(α1)-1)/np.cos(α1))*np.trapz([v_MTU1(X[:,i]) for i in range(np.shape(X)[1])],Time)
    return(error)

def muscle_length_1(X,T,T_i,lm_i,Time):
    assert np.shape(X)[0]>=2, "X must be an array of shape (M,N) where M>=2."
    lm = (1/np.cos(α1))*np.trapz([v_MTU1(X[:,i]) for i in range(np.shape(X)[1])],Time) + lm_i - (lTo1*kT/np.cos(α1))*np.log((np.exp(T/(F_MAX1*cT*kT))-1)/(np.exp(T_i/(F_MAX1*cT*kT))-1))
    return(lm)

def return_required_tension(Omega,T_i):
    Theta_i = np.pi/6
    Theta_f = 2*np.pi/3
    T_end = (Theta_f-Theta_i)/Omega
    # N = int(T_end*1000+1)
    N = 1001
    Time = np.linspace(0,T_end,N)
    X = np.zeros((2,len(Time)))
    X[0,:] = Theta_i + Omega*Time
    X[1,:] = Omega*np.ones(np.shape(Time))
    Tension = np.zeros(np.shape(Time))
    statusbar = dsb(0,len(Time),title="Finding Tensions")
    Tension = (F_MAX1*cT*kT)*np.log(np.exp((1/(lTo1*kT))*cumtrapz([v_MTU1(X[:,i]) for i in range(len(Time))],Time))*(np.exp(T_i/(F_MAX1*cT*kT))-1) + 1)
    return(Tension,Time)

def return_required_activation(T):
    u = (T*np.cos(α1) - F_MAX1*np.cos(α1)**2*F_PE1_1([0,0,0,0,lo1,0,0,0]))/(F_MAX1*np.cos(α1)**2)
    return(u)


U1 = []
T1 = []
Time1 = []
Omega1 = np.arange(0.01,0.51,0.01)
T_i1 = 100*np.ones(np.shape(Omega1))
statusbar = dsb(0,len(Omega1),title="Fixed T_i, Sweeping Omega")
for i in range(len(Omega1)):
    T_temp,Time_temp = return_required_tension(Omega1[i],T_i1[i])
    T1.append(T_temp)
    Time1.append(Time_temp)
    U1.append(return_required_activation(T_temp))
    statusbar.update(i)

plt.figure()
ax1 = plt.gca()
ax1.set_title(r"Fixed $T_{i}$, Sweeping $\omega$" + "\n Tension vs. Time")
ax1.set_xlabel("Time (s)")
ax1.set_ylabel("Tension (T)")
for i in range(len(Omega1)):
    ax1.plot(Time1[i][:-1],T1[i])

plt.figure()
ax2 = plt.gca()
ax2.set_title(r"Fixed $T_{i}$, Sweeping $\omega$" + "\n Activation vs. Time")
ax2.set_xlabel("Time (s)")
ax2.set_ylabel("Activation")
for i in range(len(Omega1)):
    ax2.plot(Time1[i][:-1],U1[i])

U2 = []
T2 = []
Time2 = []
Omega2 = 0.01*np.ones(np.shape(Omega1))
T_i2 = np.linspace(10,F_MAX1,len(Omega1))
statusbar = dsb(0,len(T_i2),title="Fixed Omega, Sweeping T_i")
for i in range(len(T_i2)):
    T_temp,Time_temp = return_required_tension(Omega2[i],T_i2[i])
    T2.append(T_temp)
    Time2.append(Time_temp)
    U2.append(return_required_activation(T_temp))
    statusbar.update(i)

plt.figure()
ax3 = plt.gca()
ax3.set_title(r"Fixed $\omega$, Sweeping $T_{i}$" + "\n Tension vs. Time")
ax3.set_xlabel("Time (s)")
ax3.set_ylabel("Tension (T)")
for i in range(len(Omega2)):
    ax3.plot(Time2[i][:-1],T2[i])

plt.figure()
ax4 = plt.gca()
ax4.set_title(r"Fixed $\omega$, Sweeping $T_{i}$" + "\n Tension vs. Time")
ax4.set_xlabel("Time (s)")
ax4.set_ylabel("Activation")
for i in range(len(Omega2)):
    ax4.plot(Time2[i][:-1],U2[i])

plt.show()

T_array = np.linspace(
    F_MAX1*np.cos(α1)*F_PE1_1([0,0,0,0,lo1,0,0,0])+0.0001,
    F_MAX1*np.cos(α1)*(1+F_PE1_1([0,0,0,0,lo1,0,0,0])),
    1001
)
plt.plot(
    T_array,
    [np.log((np.exp(T/(F_MAX1*cT*kT))-1)/(np.exp(T_i/(F_MAX1*cT*kT))-1)) for T in T_array]
)
plt.plot(Time1[0][:-1],np.exp((1/(lTo1*kT))*cumtrapz([v_MTU1(X[:,i]) for i in range(len(Time1[0]))],Time1[0])))
plt.show()
