import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

# define CSTR model
def cstr(x,t,Tc):
    Ca = x[0]
    T = x[1]
    Tf = 350
    Caf = 1.0
    q = 100
    V = 100
    rho = 1000
    Cp = 0.239
    mdelH = 5e4
    EoverR = 8750
    k0 = 7.2e10
    UA = 5e4
    rA = k0*np.exp(-EoverR/T)*Ca
    dCadt = q/V*(Caf - Ca) - rA
    dTdt = q/V*(Tf - T) \
            + mdelH/(rho*Cp)*rA \
            + UA/V/rho/Cp*(Tc-T)
    xdot = np.zeros(2)
    xdot[0] = dCadt
    xdot[1] = dTdt
    return xdot

# Steady State Initial Conditions for the States
Ca_ss = 0.87725294608097
T_ss = 324.475443431599
x0 = np.empty(2)
x0[0] = Ca_ss
x0[1] = T_ss

# Steady State Initial Condition
Tc_ss = 300.0

# Time Interval (min)
t = np.linspace(0,50,501)

# Store results for plotting
Ca = np.ones(len(t)) * Ca_ss
T = np.ones(len(t)) * T_ss
Tc = np.ones(len(t)) * Tc_ss

# Step cooling temperature
Tc[10:100]  = 303.0
Tc[100:200] = 297.0
Tc[200:300] = 300.0
Tc[300:350] = 290.0
Tc[350:400] = 302.0
Tc[400:450] = 302.0
Tc[450:]    = 299.0

# Simulate CSTR
for i in range(len(t)-1):
    ts = [t[i],t[i+1]]
    y = odeint(cstr,x0,ts,args=(Tc[i+1],))
    Ca[i+1] = y[-1][0]
    T[i+1] = y[-1][1]
    x0[0] = Ca[i+1]
    x0[1] = T[i+1]

# Construct results and save data file
# Column 1 = time
# Column 2 = cooling temperature
# Column 3 = reactor temperature
data = np.vstack((t,Tc,T)) # vertical stack
data = data.T              # transpose data
np.savetxt('cstr_step_tests.txt',data,delimiter=',',\
           header='Time,Tc,T',comments='')

# Plot the results
plt.figure()
plt.subplot(3,1,1)
plt.plot(t,Tc,'b--',linewidth=3)
plt.ylabel('Cooling T (K)')
plt.legend(['Jacket Temperature'],loc='best')

plt.subplot(3,1,2)
plt.plot(t,Ca,'r-',linewidth=3)
plt.ylabel('Ca (mol/L)')
plt.legend(['Reactor Concentration'],loc='best')

plt.subplot(3,1,3)
plt.plot(t,T,'k.-',linewidth=3)
plt.ylabel('T (K)')
plt.xlabel('Time (min)')
plt.legend(['Reactor Temperature'],loc='best')

plt.show()
