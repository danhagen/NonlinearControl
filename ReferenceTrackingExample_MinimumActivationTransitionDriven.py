import numpy as np
import matplotlib.pyplot as plt
import time
from scipy import integrate
from termcolor import cprint,colored
import matplotlib._pylab_helpers
from danpy.sb import dsb,get_terminal_width
from muscle_settings import *
from muscle_params_BIC_TRI_5_6 import *
from pendulum_equations_1DOF_2DOA import *
from pendulum_state_equations import *

"""

################################
###### Activation Driven #######
################################

x_1 &= \theta \\
x_2 &= \dot{\theta} \\
x_3 &= T_{1} \\
x_4 &= T_{2} \\
x_5 &= l_{m,1} \\
x_6 &= l_{m,2} \\
x_7 &= v_{m,1} \\
x_8 &= v_{m,2} \\
u_1 &= \alpha_1 \\
u_2 &= \alpha_2 \\

"""
N_seconds = 1

Amp = 7.5*np.pi/180
Base = 90*np.pi/180
Freq = 2*np.pi

k1,k2,k3,k4 = 100,100,100,10

if g == 0:
	MaxStep_Tension = 0.03 # percentage of positive maximum.
	Tension_Bounds = [[0,F_MAX1],[0,F_MAX2]]

	MaxStep_MuscleVelocity = 0.05 # percentage of positive maximum.
	MuscleVelocity_Bounds =[[-5*lo1,5*lo1],[-5*lo2,5*lo2]]
else:
	MaxStep_Tension = 0.01 # percentage of positive maximum.
	Tension_Bounds = [[0,F_MAX1],[0,0.10*F_MAX2]]

	MaxStep_MuscleVelocity = 1 # percentage of positive maximum.
	MuscleVelocity_Bounds =[[-5*lo1,5*lo1],[-1*lo2,1*lo2]]

MaxStep_Activation = 0.003125 # percentage of positive maximum (1)
Activation_Bounds = [[0,1],[0,1]]

"""
c_{1} &= -\frac{3g}{2L} \\
c_{2} &= \frac{3}{ML^2} \\
c_{3} &= \cos(\rho_1) \\
c_{4} &= \cos(\rho_2)
c_{5} &= \frac{\cos(\alpha_{1})}{m_1} \\
c_{6} &= \frac{\cos^2(\alpha_{1})}{m_1} \\
c_{7} &= \frac{b_{m,1}\cos^2(\alpha_{1})}{m_1} \\
c_{8} &= \tan^2(\alpha_{1}) \\
c_{9} &= \frac{\cos(\alpha_{2})}{m_2} \\
c_{10} &= \frac{\cos^2(\alpha_{2})}{m_2} \\
c_{11} &= \frac{b_{m,2}\cos^2(\alpha_{2})}{m_2} \\
c_{12} &= \tan^2(\alpha_{2}) \\

"""


# '''
# R_{1} &= r_{1}\left(\theta,\theta_{PS}=\frac{\pi}{2}\right) \\
# R_{2} &= r_{2}\left(\theta,\theta_{PS}=\frac{\pi}{2}\right)  \\
# K_{T,1} &= \frac{F_{max,1}c^{T}}{l_{T,o,1}}\left(1 - \exp{\left(\frac{-T_1}{F_{max,1}c^{T}k^{T}}\right)}\right) \hspace{1em} \text{(Variable stiffness coefficient from Zajac (1989) ODE for tendon force)} \\
# v_{MTU,1} &= \text{sgn}\left(-r_1(\theta)\right)\cdot\dot{\theta}\cdot\sqrt{\left(\frac{\partial r_1}{\partial\theta}\right)^2 + r_1^2(\theta)} \\
# K_{T,2} &= \frac{F_{max,2}c^{T}}{l_{T,o,2}}\left(1 - \exp{\left(\frac{-T_2}{F_{max,2}c^{T}k^{T}}\right)}\right) \hspace{1em} \text{(Variable stiffness coefficient from Zajac (1989) ODE for tendon force)} \\
# v_{MTU,2} &= \text{sgn}\left(-r_2(\theta)\right)\cdot\dot{\theta}\cdot\sqrt{\left(\frac{\partial r_2}{\partial\theta}\right)^2 + r_2^2(\theta)} \\
# F_{LV,1} &= f_{L,1}(l_{m,1}) \cdot f_{V,1}(l_{m,1},v_{m,1}) \\
# F_{LV,2} &= f_{L,2}(l_{m,2}) \cdot f_{V,2}(l_{m,2},v_{m,2}) \\
# '''
#
# FL = lambda l,lo: np.exp(-abs(((l/lo)**β-1)/ω)**ρ)
# FV = lambda l,v,lo: np.piecewise(v,[v<=0, v>0],\
# 	[lambda v: (V_max - v/lo)/(V_max + (cv0 + cv1*(l/lo))*(v/lo)),\
# 	lambda v: (bv-(av0 + av1*(l/lo) + av2*(l/lo)**2)*(v/lo))/(bv + (v/lo))])
#
# def R1(X):
# 	return(r1(X[0])) #
# def dR1_dx1(X):
# 	return(dr1_dθ(X[0]))
# def d2R1_dx12(X):
# 	return(d2r1_dθ2(X[0]))
# def R2(X):
# 	return(r2(X[0])) #
# def dR2_dx1(X):
# 	return(dr2_dθ(X[0]))
# def d2R2_dx12(X):
# 	return(d2r2_dθ2(X[0]))
# def KT_1(X):
# 	return((F_MAX1*cT/lTo1)*(1-np.exp(-X[2]/(F_MAX1*cT*kT)))) # NOT NORMALIZED (in N/m)
# def dKT_1_dx3(X):
# 	return((1/(kT*lTo1))*np.exp(-X[2]/(F_MAX1*cT*kT))) # NOT NORMALIZED (in N/m)
# def v_MTU1(X):
# 	return(np.sign(-R1(X))*X[1]*np.sqrt(dR1_dx1(X)**2 + R1(X)**2)) # NOT NORMALIZED (in m/s)
# def a_MTU1(X):
# 	return(np.sign(-R1(X))*(dX2_dt(X)*np.sqrt(dR1_dx1(X)**2 + R1(X)**2) \
# 				+ (X[1]**2)*dR1_dx1(X)*(d2R1_dx12(X) + R1(X))/np.sqrt(dR1_dx1(X)**2 + R1(X)**2)))
# def KT_2(X):
# 	return((F_MAX2*cT/lTo2)*(1-np.exp(-X[3]/(F_MAX2*cT*kT)))) # NOT NORMALIZED (in N/m)
# def dKT_2_dx4(X):
# 	return((1/(kT*lTo2))*np.exp(-X[3]/(F_MAX2*cT*kT))) # NOT NORMALIZED (in N/m)
# def v_MTU2(X):
# 	return(np.sign(-R2(X))*X[1]*np.sqrt(dR2_dx1(X)**2 + R2(X)**2)) # NOT NORMALIZED (in m/s)
# def a_MTU2(X):
# 	return(np.sign(-R2(X))*(dX2_dt(X)*np.sqrt(dR2_dx1(X)**2 + R2(X)**2) \
# 				+ (X[1]**2)*dR2_dx1(X)*(d2R2_dx12(X) + R2(X))/np.sqrt(dR2_dx1(X)**2 + R2(X)**2)))
# def FLV_1(X):
# 	return(FL(X[4],lo1)*FV(X[4],X[6],lo1))
# def FLV_2(X):
# 	return(FL(X[5],lo2)*FV(X[5],X[7],lo2))
# def F_PE1_1(X):
# 	return(c_1*k_1*np.log(np.exp((X[4]/(lo1*L_CE_max_1) - Lr1)/k_1) + 1) + η*(X[6]/lo1))
# def F_PE1_2(X):
# 	return(c_1*k_1*np.log(np.exp((X[5]/(lo2*L_CE_max_2) - Lr1)/k_1) + 1) + η*(X[7]/lo2))

"""
################################
######## Tension Driven ########
################################

\dot{x}_1 &= x_{2} \\
\dot{x}_2 &= c_{1}\sin(x_{1}) + c_{2}R_{1}u_{1} - c_{2}R_{2}u_{2} \\
u_1 &= T_{1} \\
u_2 &= T_{2} \\

################################
#### Muscle Velocity Driven ####
################################

\dot{x}_1 &= x_{2} \\
\dot{x}_2 &= c_{1}\sin(x_{1}) + c_{2}R_{1}x_{3} - c_{2}R_{2}x_{4} \\
\dot{x}_3 &= K_{T,1}(v_{MTU,1} - c_{3}u_1) \\
\dot{x}_4 &= K_{T,2}(v_{MTU,2} - c_{4}u_2) \\
u_1 &= \dot{l}_{m,1} \\
u_2 &= \dot{l}_{m,2} \\

################################
### Muscle Activation Driven ###
################################

\dot{x}_1 &= x_{2} \\
\dot{x}_2 &= c_{1}\sin(x_{1}) + c_{2}R_{1}x_{3} - c_{2}R_{2}x_{4} \\
\dot{x}_3 &= K_{T,1}(v_{MTU,1} - c_{3}u_1) \\
\dot{x}_4 &= K_{T,2}(v_{MTU,2} - c_{4}u_2) \\
\dot{x}_5 &= x_7 \\
\dot{x}_6 &= x_8 \\
\dot{x}_7 &= c_5x_3 - c_6F_{PE,1}(x_5,x_7) - c_7x_7 + \frac{c_{8}x_7^2}{x_5} - c_6F_{LV,1}(x_5,x_7)u_1 \\
\dot{x}_8 &= c_9x_4 - c_{10}F_{PE,2}(x_6,x_8) - c_{11}x_8 + \frac{c_{12}x_8^2}{x_6} - c_{10}F_{LV,2}(x_6,x_8)u_2 \\
u_1 &= \alpha_{1} \\
u_2 &= \alpha_{2} \\

"""
#
# def dX1_dt(X):
# 	return(X[1])
# def d2X1_dt2(X):
# 	return(dX2_dt(X))
# def dX2_dt(X,U=None):
# 	if U is None:
# 		return(c1*np.sin(X[0]) + c2*R1(X)*X[2] + c2*R2(X)*X[3])
# 	else:
# 		return(c1*np.sin(X[0]) + c2*R1(X)*U[0] + c2*R2(X)*U[1])
# def d2X2_dt2(X):
# 	return(c1*np.cos(X[0])*dX1_dt(X) + c2*dR1_dx1(X)*dX1_dt(X)*X[2] + c2*R1(X)*dX3_dt(X)\
# 			+ c2*dR2_dx1(X)*dX1_dt(X)*X[3] + c2*R2(X)*dX4_dt(X))
# def dX3_dt(X,U=None):
# 	if U is None:
# 		return(KT_1(X)*(v_MTU1(X) - c3*X[6]))
# 	else:
# 		return(KT_1(X)*(v_MTU1(X) - c3*U[0]))
# def dX4_dt(X,U=None):
# 	if U is None:
# 		return(KT_2(X)*(v_MTU2(X) - c4*X[7]))
# 	else:
# 		return(KT_2(X)*(v_MTU2(X) - c4*U[1]))
# def dX5_dt(X):
# 	return(X[6])
# def dX6_dt(X):
# 	return(X[7])
# def dX7_dt(X,U):
# 	return(c5*X[2] - c6*F_PE1_1(X) - c7*X[6] + c8*X[6]**2/X[4] - c6*FLV_1(X)*U[0])
# def dX8_dt(X,U):
# 	return(c9*X[3] - c10*F_PE1_2(X) - c11*X[7] + c12*X[7]**2/X[5] - c10*FLV_2(X)*U[1])


### Reference Trajectory ###

r = lambda t: Amp*np.cos(Freq*t) + Base
dr = lambda t: -Amp*Freq*np.sin(Freq*t)
d2r = lambda t: -Amp*Freq**2*np.cos(Freq*t)
d3r = lambda t: Amp*Freq**3*np.sin(Freq*t)
d4r = lambda t: Amp*Freq**4*np.cos(Freq*t)

############################

def Z1(t,X):
	return(r(t) - X[0])
def dZ1(t,X):
	return(dr(t) - dX1_dt(X))
def d2Z1(t,X):
	return(d2r(t) - dX2_dt(X))
def d3Z1(t,X):
	return(d3r(t) - d2X2_dt2(X))
def A1(t,X):
	return(dr(t) + k1*Z1(t,X))
def dA1(t,X):
	return(d2r(t) + k1*dZ1(t,X))
def d2A1(t,X):
	return(d3r(t) + k1*d2Z1(t,X))
def d3A1(t,X):
	return(d4r(t) + k1*d3Z1(t,X))
def Z2(t,X):
	return(X[1] - A1(t,X))
def dZ2(t,X):
	"""
	dZ2(t,X,U) = c1*np.sin(X[0]) + c2*R1(X)*U[0] + c2*R2(X)*U[1] - dA1(t,X)
	"""
	return(dX2_dt(X) - dA1(t,X))
def d2Z2(t,X):
	return(d2X2_dt2(X) - d2A1(t,X))
def A2(t,X):
	return(Z1(t,X) + dA1(t,X) - c1*np.sin(X[0]) - k2*Z2(t,X))
def dA2(t,X):
	return(dZ1(t,X) + d2A1(t,X) - c1*np.cos(X[0])*dX1_dt(X) - k2*dZ2(t,X))
def d2A2(t,X):
	return(d2Z1(t,X) + d3A1(t,X) + c1*np.sin(X[0])*(dX1_dt(X)**2) - c1*np.cos(X[0])*d2X1_dt2(X) - k2*d2Z2(t,X))
def Z3(t,X):
	return(c2*R1(X)*X[2] + c2*R2(X)*X[3] - A2(t,X))
def dZ3(t,X):
	"""
	dZ3(t,X) = c2*dR1_dx1(X)*X[1]*X[2] + c2*dR2_dx1(X)*X[1]*X[3] \
						+ c2*R1(X)*KT_1(X)*v_MTU1(X) - c2*c3*R1(X)*KT_1(X)*U[0] \
							+ c2*R2(X)*KT_2(X)*v_MTU2(X) - c2*c4*R2(X)*KT_2(X)*U[1] \
								- dA2(t,X)
	"""
	return(c2*dR1_dx1(X)*X[1]*X[2] + c2*dR2_dx1(X)*X[1]*X[3] \
					+ c2*R1(X)*KT_1(X)*v_MTU1(X) - c2*c3*R1(X)*KT_1(X)*X[6] \
						+ c2*R2(X)*KT_2(X)*v_MTU2(X) - c2*c4*R2(X)*KT_2(X)*X[7] \
							- dA2(t,X))
def A3(t,X):
	return(Z2(t,X) - dA2(t,X) + k3*Z3(t,X) \
		+ c2*dR1_dx1(X)*dX1_dt(X)*X[2] + 	c2*dR2_dx1(X)*dX1_dt(X)*X[3] \
			+ c2*R1(X)*KT_1(X)*v_MTU1(X) + c2*R2(X)*KT_2(X)*v_MTU2(X))
def dA3(t,X):
	return(dZ2(t,X) - d2A2(t,X) + k3*dZ3(t,X) \
		+ c2*d2R1_dx12(X)*(dX1_dt(X)**2)*X[2] \
			+ c2*dR1_dx1(X)*d2X1_dt2(X)*X[2] \
	 			+ c2*dR1_dx1(X)*dX1_dt(X)*dX3_dt(X)\
		 + c2*d2R2_dx12(X)*(dX1_dt(X)**2)*X[3] \
	 		+ c2*dR2_dx1(X)*d2X1_dt2(X)*X[3] \
	  			+ c2*dR2_dx1(X)*dX1_dt(X)*dX4_dt(X) \
		+ c2*dR1_dx1(X)*dX1_dt(X)*KT_1(X)*v_MTU1(X) \
			+ c2*R1(X)*dKT_1_dx3(X)*dX3_dt(X)*v_MTU1(X) \
			 	+ c2*R1(X)*KT_1(X)*a_MTU1(X) \
		+ c2*dR2_dx1(X)*dX1_dt(X)*KT_2(X)*v_MTU2(X) \
			+ c2*R2(X)*dKT_2_dx4(X)*dX4_dt(X)*v_MTU2(X) \
				+ c2*R2(X)*KT_2(X)*a_MTU2(X))
def Z4(t,X):
	return(c2*c3*R1(X)*KT_1(X)*X[6] + c2*c4*R2(X)*KT_2(X)*X[7] - A3(t,X))
def dZ4(t,X,U):
	"""
	dZ4 = 	c2*c3*dR1_dx1(X)*dX1_dt(X)*KT_1(X)*X[6]\
				+ c2*c3*R1(X)*dKT_1_dx3(X)*dX3_dt(X)*X[6]\
					+ c2*c3*R1(X)*KT_1(X)*dX7_dt(X)\
			+ c2*c4*dR2_dx1(X)*dX1_dt(X)*KT_2(X)*X[7]\
				+ c2*c4*R2(X)*dKT_2_dx4(X)*dX4_dt(X)*X[7]\
					+ c2*c4*R2(X)*KT_2(X)*dX8_dt(X)\
			- dA3(t,X)
	"""
	return(	c2*c3*dR1_dx1(X)*dX1_dt(X)*KT_1(X)*X[6] \
				+ c2*c3*R1(X)*dKT_1_dx3(X)*dX3_dt(X)*X[6] \
					+ c2*c3*R1(X)*KT_1(X)*(c5*X[2] - c6*F_PE1_1(X) - c7*X[6] + c8*X[6]**2/X[4]) \
					 	- c2*c3*c6*R1(X)*KT_1(X)*FLV_1(X)*U[0] \
			+ c2*c4*dR2_dx1(X)*dX1_dt(X)*KT_2(X)*X[7] \
			 	+ c2*c4*R2(X)*dKT_2_dx4(X)*dX4_dt(X)*X[7] \
					+ c2*c4*R2(X)*KT_2(X)*(c9*X[3] - c10*F_PE1_2(X) - c11*X[7] + c12*X[7]**2/X[5]) \
					 	- c2*c4*c10*R2(X)*KT_2(X)*FLV_2(X)*U[1] \
			- dA3(t,X))
def A4(t,X):
	"""
	c2*c3*c6*R1(X)*KT_1(X)*FLV_1(X)*U[0] \
		+ c2*c4*c10*R2(X)*KT_2(X)*FLV_2(X)*U[1] = \
					c2*c3*dR1_dx1(X)*dX1_dt(X)*KT_1(X)*X[6] \
					+ c2*c3*R1(X)*dKT_1_dx3(X)*dX3_dt(X)*X[6] \
					+ c2*c3*R1(X)*KT_1(X)*(c5*X[2] - c6*F_PE1_1(X) - c7*X[6] + c8*X[6]**2/X[4]) \
					+ c2*c4*dR2_dx1(X)*dX1_dt(X)*KT_2(X)*X[7] \
					+ c2*c4*R2(X)*dKT_2_dx4(X)*dX4_dt(X)*X[7] \
					+ c2*c4*R2(X)*KT_2(X)*(c9*X[3] - c10*F_PE1_2(X) - c11*X[7] + c12*X[7]**2/X[5]) \
					- dA3(t,X) - Z3(t,X) + k4*Z4(t,X)
	"""
	return(c2*c3*dR1_dx1(X)*dX1_dt(X)*KT_1(X)*X[6] \
				+ c2*c3*R1(X)*dKT_1_dx3(X)*dX3_dt(X)*X[6] \
					+ c2*c3*R1(X)*KT_1(X)*(c5*X[2] - c6*F_PE1_1(X) - c7*X[6] + c8*X[6]**2/X[4]) \
			+ c2*c4*dR2_dx1(X)*dX1_dt(X)*KT_2(X)*X[7] \
			 	+ c2*c4*R2(X)*dKT_2_dx4(X)*dX4_dt(X)*X[7] \
					+ c2*c4*R2(X)*KT_2(X)*(c9*X[3] - c10*F_PE1_2(X) - c11*X[7] + c12*X[7]**2/X[5]) \
			- dA3(t,X) - Z3(t,X) + k4*Z4(t,X))

def return_constraint_variables_muscle_activation_driven(t,X):
	Coefficient1 = c2*c3*c6*R1(X)*KT_1(X)*FLV_1(X)
	Coefficient2 = c2*c4*c10*R2(X)*KT_2(X)*FLV_2(X)
	Constraint = A4(t,X)
	return(Coefficient1,Coefficient2,Constraint)

def return_initial_tension(X_o,**kwargs):
	"""
	Takes in initial state numpy.ndarray (X_o) of shape (2,) and returns an initial tension (2,) numpy.ndarray for .

	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	**kwargs
	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	1) Bounds - must be a (2,2) list with each row in ascending order. Default is given by Tension_Bounds.

	2) InitialAngularAcceleration - must be a float or an int. Default is 0 (starting from rest).

	"""

	assert (np.shape(X_o) in [(2,),(4,),(8,)]) \
			and (str(type(X_o)) == "<class 'numpy.ndarray'>"), \
		"X_o must be a (2,), (4,), or (8,) numpy.ndarray."

	InitialAngularAcceleration = kwargs.get("InitialAngularAcceleration",0) # or d2r(0)
	assert type(InitialAngularAcceleration) in [float,int], "InitialAngularAcceleration must be either a float or an int."

	Bounds = kwargs.get("Bounds",Tension_Bounds)
	assert type(Bounds) == list and np.shape(Bounds) == (2,2), "Bounds for Tension Control must be a (2,2) list."
	assert Bounds[0][0]<Bounds[0][1],"Each set of bounds must be in ascending order."
	assert Bounds[1][0]<Bounds[1][1],"Each set of bounds must be in ascending order."

	Constraint = lambda T1: -(R1(X_o)*T1 + (M*L**2/3)*InitialAngularAcceleration)/R2(X_o)
	InverseConstraint = lambda T2: -(R2(X_o)*T2 + (M*L**2/3)*InitialAngularAcceleration)/R1(X_o)

	LowerBound_x = max(Bounds[0][0],InverseConstraint(Bounds[1][0]))
	LowerBound_y = Constraint(LowerBound_x)
	UpperBound_x = min(Bounds[0][1],InverseConstraint(Bounds[1][1]))
	UpperBound_y = Constraint(UpperBound_x)

	LowerBoundVector = np.array([[LowerBound_x],[LowerBound_y]])
	UpperBoundVector = np.array([[UpperBound_x],[UpperBound_y]])

	InitialTension = (UpperBoundVector-LowerBoundVector)*np.random.rand() + LowerBoundVector
	return(InitialTension)
def plot_initial_tension_values(X_o,**kwargs):
	"""
	Takes in initial state numpy.ndarray (X_o) of shape (2,) and plots 1000 initial tension values (2,).

	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	**kwargs
	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	1) Seed - must be a scalar value. Default is None.

	2) Bounds - must be a (2,2) list with each row in ascending order. Default is given by Tension_Bounds.

	3) InitialAngularAcceleration - must be a float or an int. Default is 0 (starting from rest).

	4) NumberOfPoints - must be an int. Default is 1000.

	"""

	assert (np.shape(X_o) in [(2,),(4,),(8,)]) \
			and (str(type(X_o)) == "<class 'numpy.ndarray'>"), \
		"X_o must be a (2,), (4,), or (8,) numpy.ndarray."

	InitialAngularAcceleration = kwargs.get("InitialAngularAcceleration",0) # or d2r(0)
	assert type(InitialAngularAcceleration) in [float,int], "InitialAngularAcceleration must be either a float or an int."

	Seed = kwargs.get("Seed",None)
	assert type(Seed) in [float,int] or Seed is None, "Seed must be a float or an int or None."
	np.random.seed(Seed)

	Bounds = kwargs.get("Bounds",Tension_Bounds)
	assert type(Bounds) == list and np.shape(Bounds) == (2,2), "Bounds for Tension Control must be a (2,2) list."
	assert Bounds[0][0]<Bounds[0][1],"Each set of bounds must be in ascending order."
	assert Bounds[1][0]<Bounds[1][1],"Each set of bounds must be in ascending order."

	NumberOfPoints = kwargs.get("NumberOfPoints",1000)
	assert type(NumberOfPoints) == int, "NumberOfPoints must be an int."

	fig = plt.figure(figsize=(10,8))
	ax1 = plt.gca()

	DescriptiveTitle = "Plotting Randomly Generated\nInitial Tendon Tensions"

	ax1.set_title(DescriptiveTitle,Fontsize=20,y=1)

	#Bound Constraints
	ax1.plot([Bounds[0][0],Bounds[0][1]],[Bounds[1][0],Bounds[1][0]],'k--')
	ax1.plot([Bounds[0][0],Bounds[0][1]],[Bounds[1][1],Bounds[1][1]],'k--')
	ax1.plot([Bounds[0][0],Bounds[0][0]],[Bounds[1][0],Bounds[1][1]],'k--')
	ax1.plot([Bounds[0][1],Bounds[0][1]],[Bounds[1][0],Bounds[1][1]],'k--')

	Constraint = lambda T1: -(R1(X_o)*T1 + (M*L**2/3)*InitialAngularAcceleration)/R2(X_o)
	Input1 = np.linspace(
		Bounds[0][0]-0.2*np.diff(Bounds[0]),
		Bounds[0][1]+0.2*np.diff(Bounds[0]),
		1001)
	Input2 = Constraint(Input1)
	ax1.plot(Input1,Input2,'r',lw=2)
	statusbar = dsb()
	for i in range(NumberOfPoints):
		T = return_initial_tension(X_o,**kwargs)
		ax1.plot(T[0],T[1],'go')
		statusbar.update(i,NumberOfPoints,title=plot_initial_tension_values.__name__)
	ax1.set_xlabel(r'$T_{1}$',fontsize=14)
	ax1.set_ylabel(r'$T_{2}$',fontsize=14)
	ax1.set_xlim([Input1[0],Input1[-1]])
	ax1.set_ylim([Input2[0],Input2[-1]])
	ax1.spines['right'].set_visible(False)
	ax1.spines['top'].set_visible(False)
	ax1.set_aspect('equal')
	plt.show()
def return_initial_muscle_lengths_and_activations(InitialTension,**kwargs):
	Statusbar = kwargs.get("Statusbar",False)
	assert type(Statusbar)==bool, "Statusbar must be a boolean. Default is False."

	PlotBool = kwargs.get("PlotBool",False)
	assert type(PlotBool)==bool,"PlotBool must be a boolean. Default is False."

	# L1 = np.linspace(0.5*lo1, 1.5*lo1, 1001)
	mu1, sigma1 = lo1, 0.1*lo1
	L1 = np.array(list(sorted(np.random.normal(mu1, sigma1, 1001))))
	U1 = (
		InitialTension[0][0]/(F_MAX1*np.cos(α1))
	    - c_1*k_1*np.log(np.exp((L1/(lo1*L_CE_max_1) - Lr1)/k_1)+1)
	    ) / (np.exp(-(abs((L1-lo1)/(lo1*ω))**ρ)))

	# L2 = np.linspace(0.5*lo2, 1.5*lo2, 1001)
	mu2, sigma2 = lo2, 0.1*lo2
	L2 = np.array(list(sorted(np.random.normal(mu2, sigma2, 1001))))
	U2 = (
		InitialTension[1][0]/(F_MAX2*np.cos(α2))
	    - c_1*k_1*np.log(np.exp((L2/(lo2*L_CE_max_2) - Lr1)/k_1)+1)
	    ) / (np.exp(-(abs((L2-lo2)/(lo2*ω))**ρ)))

	if PlotBool == True:
		plt.figure(figsize=(10,8))
		plt.title(r"Viable Initial $l_{m,1}$ and $u_{1}$ Values")
		plt.xlabel(r"$l_{m,1}$ (m)",fontsize=14)
		plt.ylabel(r"$u_{1}$",fontsize=14)
		plt.scatter(L1,U1)
		plt.plot([lo1,lo1],[0,1],'0.70',linestyle='--')
		plt.gca().set_ylim((0,1))
		plt.gca().set_xticks(
			[0.25*lo1,
			0.5*lo1,
			0.75*lo1,
			lo1,
			1.25*lo1,
			1.5*lo1,
			1.75*lo1]
			)
		plt.gca().set_xticklabels(
			["",
			r"$\frac{1}{2}$ $l_{o,2}$",
			"",
			r"$l_{o,2}$",
			"",
			r"$\frac{3}{2}$ $l_{o,2}$",
			""],
			fontsize=12)

		plt.figure(figsize=(10,8))
		plt.title(r"Viable Initial $l_{m,2}$ and $u_{2}$ Values")
		plt.xlabel(r"$l_{m,2}$ (m)",fontsize=14)
		plt.ylabel(r"$u_{2}$",fontsize=14)
		plt.scatter(L2,U2)
		plt.plot([lo2,lo2],[0,1],'0.70',linestyle='--')
		plt.gca().set_ylim((0,1))
		plt.gca().set_xticks(
			[0.25*lo2,
			0.5*lo2,
			0.75*lo2,
			lo2,
			1.25*lo2,
			1.5*lo2,
			1.75*lo2]
			)
		plt.gca().set_xticklabels(
			["",
			r"$\frac{1}{2}$ $l_{o,2}$",
			"",
			r"$l_{o,2}$",
			"",
			r"$\frac{3}{2}$ $l_{o,2}$",
			""],
			fontsize=12)

		plt.show()
	return(L1,U1,L2,U2)
def find_viable_initial_values(**kwargs):
	ReturnAll = kwargs.get("ReturnAll",False)
	assert type(ReturnAll)==bool, "ReturnAll must be a bool."

	Seed = kwargs.get("Seed",None)
	assert type(Seed) in [float,int] or Seed is None, "Seed must be a float or an int or None."
	np.random.seed(Seed)

	X_o = np.array([Amp+Base,0])
	T = return_initial_tension(X_o)
	L1,U1,L2,U2 = return_initial_muscle_lengths_and_activations(T,**kwargs)
	rand_index = np.random.choice(len(L1),2)

	if ReturnAll == False:
		return(
			T,
			np.array([L1[rand_index[0]],L2[rand_index[1]]]),
			np.array([U1[rand_index[0]],U2[rand_index[1]]])
			)
	else:
		return(T,L1,L2,U1,U2)
def return_U_muscle_activation_driven(i,t,X,U,**kwargs):
	"""
	Takes in current step (i), numpy.ndarray of time (t) of shape (N,), state numpy.ndarray (X) of shape (8,), and previous input numpy.ndarray (U) of shape (2,) and returns the input for this time step.

	First attempt will see what happens when the activations are restricted to the positive real domain.

	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	**kwargs
	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	1) Noise - must be an numpy.ndarray of shape (2,). Default is np.zeros((1,2)).

	2) Seed - must be a scalar value. Default is None.

	3) Bounds - must be a (2,2) list with each row in ascending order. Default is given by Activation_Bounds.

	4) MaxStep - must be a scalar (int or float). Default is MaxStep_Activation.

	"""
	import random
	import numpy as np
	assert (np.shape(t) == (len(t),)) and (str(type(t)) == "<class 'numpy.ndarray'>"),\
	 	"t must be a numpy.ndarray of shape (len(t),)."
	assert np.shape(X) == (8,) and str(type(X)) == "<class 'numpy.ndarray'>", "X must be a (8,) numpy.ndarray"
	assert np.shape(U) == (2,) and str(type(U)) == "<class 'numpy.ndarray'>", "U must be a (2,) numpy.ndarray"

	dt = t[1]-t[0]

	Noise = kwargs.get("Noise",np.zeros((2,)))
	assert np.shape(Noise) == (2,) and str(type(Noise)) == "<class 'numpy.ndarray'>", "Noise must be a (2,) numpy.ndarray"

	Seed = kwargs.get("Seed",None)
	assert type(Seed) in [float,int] or Seed is None, "Seed must be a float or an int or None."
	np.random.seed(Seed)

	Bounds = kwargs.get("Bounds",Activation_Bounds)
	assert type(Bounds) == list and np.shape(Bounds) == (2,2), "Bounds for Muscle Activation Control must be a (2,2) list."
	assert Bounds[0][0]<Bounds[0][1],"Each set of bounds must be in ascending order."
	assert Bounds[1][0]<Bounds[1][1],"Each set of bounds must be in ascending order."

	Coefficient1,Coefficient2,Constraint1 =\
	 			return_constraint_variables_muscle_activation_driven(t[i],X)
	assert Coefficient1!=0 and Coefficient2!=0, "Error with Coefficients. Shouldn't both be zero."
	if Constraint1 < 0:
		assert not(Coefficient1 > 0 and Coefficient2 > 0), "Infeasible activations. (Constraint1 < 0, Coefficient1 > 0, Coefficient2 > 0)"
	if Constraint1 > 0:
		assert not(Coefficient1 < 0 and Coefficient2 < 0), "Infeasible activations. (Constraint1 > 0, Coefficient1 < 0, Coefficient2 < 0)"

	if Coefficient1 == 0:
		LowerBound_x = max(Bounds[0][0],AllowbaleBounds_x[0])
		UpperBound_x = min(Bounds[0][1],AllowbaleBounds_x[1])
		FeasibleInput1 = (UpperBound_x-LowerBound_x)*np.random.rand(1000) + LowerBound_x
		FeasibleInput2 = np.array([Constraint1/Coefficient2]*1000)
	elif Coefficient2 == 0:
		LowerBound_y = max(Bounds[1][0],AllowableBounds_y[0])
		UpperBound_y = min(Bounds[1][1],AllowableBounds_y[1])
		FeasibleInput1 = np.array([Constraint1/Coefficient1]*1000)
		FeasibleInput2 = (UpperBound_y-LowerBound_y)*np.random.rand(1000) + LowerBound_y
	else:
		SortedBounds = np.sort([(Constraint1-Coefficient2*Bounds[1][0])/Coefficient1,\
									(Constraint1-Coefficient2*Bounds[1][1])/Coefficient1])
		LowerBound_x = max(	Bounds[0][0],\
		 					SortedBounds[0]\
						)
		UpperBound_x = min(	Bounds[0][1],\
		 					SortedBounds[1]\
						)
		# if UpperBound_x < LowerBound_x: import ipdb; ipdb.set_trace()
		assert UpperBound_x >= LowerBound_x, "Error generating bounds. Not feasible!"
		# FeasibleInput1 = (UpperBound_x-LowerBound_x)*np.random.rand(1000) + LowerBound_x
		# FeasibleInput2 = np.array([Constraint1/Coefficient2 - (Coefficient1/Coefficient2)*el \
		# 						for el in FeasibleInput1])
		FeasibleInput1 = np.ones(1000)*(Coefficient2*(Coefficient2*U[0]-Coefficient1*U[1])+Coefficient1*Constraint1)/(Coefficient1**2 + Coefficient2**2)
		FeasibleInput2 = np.ones(1000)*(Coefficient1*(-Coefficient2*U[0]+Coefficient1*U[1])+Coefficient2*Constraint1)/(Coefficient1**2 + Coefficient2**2)
		assert LowerBound_x<=FeasibleInput1[0]<=UpperBound_x, "No Feasible transition. Closest transition greater than maximum allowable transition."
	"""
	Checking to see which inputs have the appropriate allowable step size.
	"""
	# euclid_dist = np.array(list(map(lambda x,y: np.sqrt((U[0]-x)**2+(U[1]-y)**2),\
	# 								FeasibleInput1,FeasibleInput2)))
	# next_index, = np.where(euclid_dist==min(euclid_dist))
	# u1 = FeasibleInput1[next_index[0]]
	# u2 = FeasibleInput2[next_index[0]]
	u1 = FeasibleInput1[0]
	u2 = FeasibleInput2[0]
	return(np.array([u1,u2]))

def return_length_of_nonzero_array(X):
	"""
	Takes in a numpy.ndarray X of shape (m,n) and returns the length of the array that removes any trailing zeros.
	"""
	assert str(type(X))=="<class 'numpy.ndarray'>", "X should be a numpy array"
	assert np.shape(X)[1]!=1, "X should be a wide rectangular array. (m,1) is a column, therefore a nonzero X of this shape will return 1 (trivial solution). Transpose X to properly identify nonzero array length."
	assert np.shape(X)!=(1,1), "Check input. Should not be of shape (1,1) (trivial solution)."
	if (X[:,1:]!=np.zeros(np.shape(X[:,1:]))).all():
		return(np.shape(X)[1])
	else:
		return(np.argmax((X[:,1:] == np.zeros(np.shape(X[:,1:]))).sum(axis=0) == np.shape(X[:,1:])[0])+1)

def run_sim_MAT(**kwargs):
	"""
	Runs one simulation for MINIMUM ACTIVATION TRANSITION control.

	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	**kwargs
	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	1) Bounds - must be a (2,2) list with each row in ascending order. Default is given by Tension_Bounds.

	2) InitialAngularAcceleration - must be a float or an int. Default is 0 (starting from rest).

	3) thresh - must be an int. Default is 25.

	"""
	thresh = kwargs.get("thresh",25)
	assert type(thresh)==int, "thresh should be an int as it is the number of attempts the program should run before stopping."

	N = N_seconds*10000 + 1
	Time = np.linspace(0,N_seconds,N)
	dt = Time[1]-Time[0]

	AnotherIteration = True
	AttemptNumber = 1

	while AnotherIteration == True:
		X = np.zeros((8,N))
		InitialTension,InitialMuscleLengths,InitialActivations = \
			find_viable_initial_values(**kwargs)
		X[:,0] = [
			r(0),
			dr(0),
			InitialTension[0][0],
			InitialTension[1][0],
			InitialMuscleLengths[0],
			InitialMuscleLengths[1],
			0,
			0]
		U = np.zeros((2,N))
		U[:,0] = InitialActivations

		AddNoise = False
		if AddNoise == True:
		    np.random.seed(seed=None)
		    NoiseArray = np.random.normal(loc=0.0,scale=0.2,size=(2,N))
		else:
		    NoiseArray = np.zeros((2,N))

		try:
			cprint("Attempt #" + str(int(AttemptNumber)) + ":\n", 'green')
			statusbar = dsb(0,N-1,title=run_sim_MAT.__name__)
			for i in range(N-1):
				U[:,i+1] = return_U_muscle_activation_driven(i,Time,X[:,i],U[:,i],Noise = NoiseArray[:,i])
				X[:,i+1] = X[:,i] + dt*np.array([	dX1_dt(X[:,i]),\
													dX2_dt(X[:,i]),\
													dX3_dt(X[:,i]),\
													dX4_dt(X[:,i]),\
													dX5_dt(X[:,i]),\
													dX6_dt(X[:,i]),\
													dX7_dt(X[:,i],U=U[:,i+1]),\
													dX8_dt(X[:,i],U=U[:,i+1])])
				statusbar.update(i)
			AnotherIteration = False
			return(X,U)
		except:
			print('\n')
			print(" "*(get_terminal_width()\
						- len("...Attempt #" + str(int(AttemptNumber)) + " Failed. "))\
						+ colored("...Attempt #" + str(int(AttemptNumber)) + " Failed. \n",'red'))
			AttemptNumber += 1
			if AttemptNumber > thresh:
				AnotherIteration=False
				return(np.zeros((8,N)),np.zeros((2,N)))
def run_N_sim_MAT(**kwargs):
	NumberOfTrials = kwargs.get("NumberOfTrials",10)

	N = N_seconds*10000 + 1
	Time = np.linspace(0,N_seconds,N)
	dt = Time[1]-Time[0]

	TotalX = np.zeros((NumberOfTrials,8,N))
	TotalU = np.zeros((NumberOfTrials,2,N))
	TerminalWidth = get_terminal_width()

	print("\n")
	for j in range(NumberOfTrials):
		TrialTitle = "          Trial #" + str(j+1)+ "          \n"
		print(
			" "*int(TerminalWidth/2 - len(TrialTitle)/2)
			+ colored(TrialTitle,'white',attrs=["underline","bold"])
			)
		TotalX[j],TotalU[j] = run_sim_MAT(**kwargs)

	i=0
	NumberOfSuccessfulTrials = NumberOfTrials
	while i < NumberOfSuccessfulTrials:
		if (TotalX[i]==np.zeros((8,np.shape(TotalX)[2]))).all():
			TotalX = np.delete(TotalX,i,0)
			TotalU = np.delete(TotalU,i,0)
			NumberOfSuccessfulTrials-=1
			if NumberOfSuccessfulTrials==0: raise ValueError("No Successful Trials!")
		else:
			i+=1

	print(
		"Number of Desired Runs: "
		+ str(NumberOfTrials)
		+ "\n"
		+ "Number of Successful Runs: "
		+ str(NumberOfSuccessfulTrials)
		+ "\n"
	)
	return(TotalX,TotalU)

def plot_MA_values(t,X,**kwargs):
	"""
	Take the numpy.ndarray time array (t) of size (N,) and the state space numpy.ndarray (X) of size (2,N), (4,N), or (8,N), and plots the moment are values of the two muscles versus time and along the moment arm function.

	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	**kwargs
	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	1) InputString - must be a string. Used to alter the figure Title. Default is None.
	"""
	import matplotlib.pyplot as plt
	import numpy as np

	assert (np.shape(X)[0] in [2,4,8]) \
				and (np.shape(X)[1] == len(t)) \
					and (str(type(X)) == "<class 'numpy.ndarray'>"), \
			"X must be a (2,N), (4,N), or (8,N) numpy.ndarray, where N is the length of t."

	assert np.shape(t) == (len(t),) and str(type(t)) == "<class 'numpy.ndarray'>", "t must be a (N,) numpy.ndarray."

	InputString = kwargs.get("InputString",None)
	assert InputString is None or type(InputString)==str, "InputString must either be a string or None."
	if InputString is None:
		DescriptiveTitle = "Moment arm equations"
	else:
		assert type(InputString)==str, "InputString must be a string"
		DescriptiveTitle = "Moment arm equations\n(" + InputString + " Driven)"

	fig, ((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2,figsize=(8,6))
	plt.subplots_adjust(left = 0.15,hspace=0.1,bottom=0.1)
	plt.suptitle(DescriptiveTitle)

	ax1.plot(np.linspace(0,np.pi*(160/180),1001),\
				np.array(list(map(lambda x1: R1([x1]),np.linspace(0,np.pi*(160/180),1001)))),\
				'0.70')
	ax1.plot(np.linspace(min(X[0,:]),max(X[0,:]),101),\
				np.array(list(map(lambda x1: R1([x1]),np.linspace(min(X[0,:]),max(X[0,:]),101)))),\
				'g',lw=3)
	ax1.set_xticks([0,np.pi/4,np.pi/2,3*np.pi/4,np.pi])
	ax1.set_xticklabels([""]*len(ax1.get_xticks()))
	ax1.set_ylabel("Moment Arm for\n Muscle 1 (m)")

	"""
	Note: Need to Transpose X in order for Map to work.
	"""

	ax2.plot(t,np.array(list(map(lambda X: R1(X),X.T))),'g')
	ax2.set_ylim(ax1.get_ylim())
	ax2.set_yticks(ax1.get_yticks())
	ax2.set_yticklabels([""]*len(ax1.get_yticks()))
	ax2.set_xticklabels([""]*len(ax2.get_xticks()))

	ax3.plot(np.linspace(0,np.pi*(160/180),1001),\
				np.array(list(map(lambda x1: R2([x1]),np.linspace(0,np.pi*(160/180),1001)))),\
				'0.70')
	ax3.plot(np.linspace(min(X[0,:]),max(X[0,:]),101),\
				np.array(list(map(lambda x1: R2([x1]),np.linspace(min(X[0,:]),max(X[0,:]),101)))),\
				'r',lw=3)
	ax3.set_xticks([0,np.pi/4,np.pi/2,3*np.pi/4,np.pi])
	ax3.set_xticklabels([r"$0$",r"$\frac{\pi}{4}$",r"$\frac{\pi}{2}$",r"$\frac{3\pi}{4}$",r"$\pi$"])
	ax3.set_xlabel("Joint Angle (rads)")
	ax3.set_ylabel("Moment Arm for\n Muscle 2 (m)")

	ax4.plot(t,np.array(list(map(lambda X: R2(X),X.T))),'r')
	ax4.set_ylim(ax3.get_ylim())
	ax4.set_yticks(ax3.get_yticks())
	ax4.set_yticklabels([""]*len(ax3.get_yticks()))
	ax4.set_xlabel("Time (s)")
	return(fig,[ax1,ax2,ax3,ax4])

def animate_muscle_activation_driven(t,X,U,**kwargs):
	"""
	Takes in Time (t - numpy.ndarray of shape (N,)), the state array (X - numpy.ndarray of shape (8,N)), and the input array (U - numpy.ndarray of shape (2,N)) and animates constraint equation over time.

	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	**kwargs
	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	1) Bounds - must be a (2,2) list with each row in ascending order. Default is given by MuscleVelocity_Bounds.

	"""
	import numpy as np
	import matplotlib.pyplot as plt
	import matplotlib.animation as animation
	import matplotlib.patches as patches
	import time

	assert np.shape(X) == (8,len(t)) and str(type(X)) == "<class 'numpy.ndarray'>", "X must be a (8,N) numpy.ndarray"

	assert np.shape(U) == (2,len(t)) and str(type(X)) == "<class 'numpy.ndarray'>", "X must be a (8,N) numpy.ndarray"

	MaxStep = kwargs.get("MaxStep",MaxStep_Activation)
	assert type(MaxStep) in [int,float], "MaxStep for Muscle Activation Controller should be an int or float."

	Bounds = kwargs.get("Bounds",Activation_Bounds)
	assert type(Bounds) == list and np.shape(Bounds) == (2,2), "Bounds for Muscle Activation Control must be a (2,2) list."
	assert Bounds[0][0]<Bounds[0][1],"Each set of bounds must be in ascending order."
	assert Bounds[1][0]<Bounds[1][1],"Each set of bounds must be in ascending order."

	dt = t[1]-t[0]
	fig = plt.figure(figsize=(10,8))
	ax1 = plt.gca()

	DescriptiveTitle = "Plotting Constraints vs. Time\nMuscle Activation Driven"

	ax1.set_title(DescriptiveTitle,Fontsize=20,y=1)

	#Bound Constraints
	ax1.plot([Bounds[0][0],Bounds[0][1]],[Bounds[1][0],Bounds[1][0]],'k--')
	ax1.plot([Bounds[0][0],Bounds[0][1]],[Bounds[1][1],Bounds[1][1]],'k--')
	ax1.plot([Bounds[0][0],Bounds[0][0]],[Bounds[1][0],Bounds[1][1]],'k--')
	ax1.plot([Bounds[0][1],Bounds[0][1]],[Bounds[1][0],Bounds[1][1]],'k--')

	Coefficient1,Coefficient2,Constraint1 = return_constraint_variables_muscle_activation_driven(t[0],X[:,0])

	AllowableBounds_x = np.array([U[0,0]-MaxStep,U[0,0]+MaxStep])
	AllowableBounds_y = np.array([U[1,0]-MaxStep,U[1,0]+MaxStep])

	if Coefficient1 == 0:
		LowerBound_x = max(Bounds[0][0],AllowbaleBounds_x[0])
		UpperBound_x = min(Bounds[0][1],AllowbaleBounds_x[1])
		LowerBound_y = Constraint1/Coefficient2
		UpperBound_y = Constraint1/Coefficient2
		FeasibleInput1 = (UpperBound_x-LowerBound_x)*np.random.rand(1000) + LowerBound_x
		FeasibleInput2 = np.array([Constraint1/Coefficient2]*1000)
	elif Coefficient2 == 0:
		LowerBound_x = Constraint1/Coefficient1
		UpperBound_x = Constraint1/Coefficient1
		LowerBound_y = max(Bounds[1][0],AllowableBounds_y[0])
		UpperBound_y = min(Bounds[1][1],AllowableBounds_y[1])
		FeasibleInput1 = np.array([Constraint1/Coefficient1]*1000)
		FeasibleInput2 = (UpperBound_y-LowerBound_y)*np.random.rand(1000) + LowerBound_y
	else:
		SortedAllowableBounds = np.sort([\
									(Constraint1-Coefficient2*AllowableBounds_y[0])/Coefficient1,\
									(Constraint1-Coefficient2*AllowableBounds_y[1])/Coefficient1\
									])
		SortedBounds = np.sort([(Constraint1-Coefficient2*Bounds[1][0])/Coefficient1,\
									(Constraint1-Coefficient2*Bounds[1][1])/Coefficient1])
		LowerBound_x = max(	Bounds[0][0],\
		 					SortedBounds[0],\
							AllowableBounds_x[0],\
							SortedAllowableBounds[0]\
						)
		UpperBound_x = min(	Bounds[0][1],\
		 					SortedBounds[1],\
							AllowableBounds_x[1],\
							SortedAllowableBounds[1]\
						)
		LowerBound_y,UpperBound_y = np.sort([\
									(Constraint1-Coefficient1*LowerBound_x)/Coefficient2,\
									(Constraint1-Coefficient1*UpperBound_x)/Coefficient2\
									])
		# if UpperBound < LowerBound: import ipdb; ipdb.set_trace()
		assert UpperBound_x >= LowerBound_x, "Error generating bounds. Not feasible!"
		FeasibleInput1 = (UpperBound_x-LowerBound_x)*np.random.rand(1000) + LowerBound_x
		FeasibleInput2 = np.array([Constraint1/Coefficient2 - (Coefficient1/Coefficient2)*el \
								for el in FeasibleInput1])

	feasible = plt.Rectangle((LowerBound_x,LowerBound_y),\
								UpperBound_x-LowerBound_x,\
								UpperBound_y-LowerBound_y,\
								alpha = 0.5,
								Color = 'b')
	# feasible = plt.Circle((U[:,0]),radius=MaxStep_Activation,Color='b',alpha=0.5)
	ax1.add_patch(feasible)
	if Coefficient2!=0:
		cline_dashed, = plt.plot(\
			ax1.get_xlim(),\
			np.array(list(map(lambda x: \
				(Constraint1-Coefficient1*x)/Coefficient2, ax1.get_xlim()))),\
			'k--',lw=1)
	else:
		cline_dashed, = plt.plot(\
			np.array(list(map(lambda y: \
				(Constraint1-Coefficient2*y)/Coefficient1, ax1.get_ylim()))),\
			ax1.get_ylim(),\
			'k--',lw=1)
	cline, = plt.plot(FeasibleInput1,FeasibleInput2,'b',lw=2)
	TimeText = plt.text(0.1,0.1,"t = " + str(t[0]),fontsize=16)
	chosenpoint, = plt.plot(U[:,0],c='k',marker='o')
	ax1.set_xlabel(r'$\alpha_{1}$',fontsize=14)
	ax1.set_ylabel(r'$\alpha_{2}$',fontsize=14)
	ax1.set_xlim([Bounds[0][0]-0.10*(np.diff(Bounds[0])[0]/2),\
					Bounds[0][1]+0.10*(np.diff(Bounds[0])[0]/2)])
	ax1.set_ylim([Bounds[1][0]-0.10*(np.diff(Bounds[1])[0]/2),\
					Bounds[1][1]+0.10*(np.diff(Bounds[1])[0]/2)])
	ax1.spines['right'].set_visible(False)
	ax1.spines['top'].set_visible(False)
	ax1.set_aspect('equal')

	def animate(i):
		Coefficient1,Coefficient2,Constraint1 = return_constraint_variables_muscle_activation_driven(t[i],X[:,i])
		AllowableBounds_x = np.array([U[0,i]-MaxStep,U[0,i]+MaxStep])
		AllowableBounds_y = np.array([U[1,i]-MaxStep,U[1,i]+MaxStep])

		if Coefficient1 == 0:
			LowerBound_x = max(Bounds[0][0],AllowbaleBounds_x[0])
			UpperBound_x = min(Bounds[0][1],AllowbaleBounds_x[1])
			LowerBound_y = Constraint1/Coefficient2
			UpperBound_y = Constraint1/Coefficient2
			FeasibleInput1 = (UpperBound_x-LowerBound_x)*np.random.rand(1000) + LowerBound_x
			FeasibleInput2 = np.array([Constraint1/Coefficient2]*1000)
		elif Coefficient2 == 0:
			LowerBound_x = Constraint1/Coefficient1
			UpperBound_x = Constraint1/Coefficient1
			LowerBound_y = max(Bounds[1][0],AllowableBounds_y[0])
			UpperBound_y = min(Bounds[1][1],AllowableBounds_y[1])
			FeasibleInput1 = np.array([Constraint1/Coefficient1]*1000)
			FeasibleInput2 = (UpperBound_y-LowerBound_y)*np.random.rand(1000) + LowerBound_y
		else:
			SortedAllowableBounds = np.sort([\
										(Constraint1-Coefficient2*AllowableBounds_y[0])/Coefficient1,\
										(Constraint1-Coefficient2*AllowableBounds_y[1])/Coefficient1\
										])
			SortedBounds = np.sort([(Constraint1-Coefficient2*Bounds[1][0])/Coefficient1,\
										(Constraint1-Coefficient2*Bounds[1][1])/Coefficient1])
			LowerBound_x = max(	Bounds[0][0],\
			 					SortedBounds[0],\
								AllowableBounds_x[0],\
								SortedAllowableBounds[0]\
							)
			UpperBound_x = min(	Bounds[0][1],\
			 					SortedBounds[1],\
								AllowableBounds_x[1],\
								SortedAllowableBounds[1]\
							)
			LowerBound_y,UpperBound_y = np.sort([\
										(Constraint1-Coefficient1*LowerBound_x)/Coefficient2,\
										(Constraint1-Coefficient1*UpperBound_x)/Coefficient2\
										])
			# if UpperBound < LowerBound: import ipdb; ipdb.set_trace()
			assert UpperBound_x >= LowerBound_x, "Error generating bounds. Not feasible!"
			FeasibleInput1 = (UpperBound_x-LowerBound_x)*np.random.rand(1000) + LowerBound_x
			FeasibleInput2 = np.array([Constraint1/Coefficient2 - (Coefficient1/Coefficient2)*el \
									for el in FeasibleInput1])
		# feasible.center = (U[:,i])
		feasible.set_xy((LowerBound_x,LowerBound_y))
		feasible.set_width(UpperBound_x-LowerBound_x)
		feasible.set_height(UpperBound_y-LowerBound_y)
		# if i<10:
		# 	feasible.radius = 10*MaxStep_Activation
		# else:
		# 	feasible.radius = MaxStep_Activation
		if Coefficient2!=0:
			cline_dashed.set_xdata(ax1.get_xlim())
			cline_dashed.set_ydata(np.array(list(map(lambda x: \
				(Constraint1-Coefficient1*x)/Coefficient2, ax1.get_xlim()))))
		else:
			cline_dashed.set_xdata(np.array(list(map(lambda y: \
				(Constraint1-Coefficient2*y)/Coefficient1, ax1.get_ylim()))))
			cline_dashed.set_ydata(ax1.get_ylim())
		cline.set_xdata(FeasibleInput1)
		cline.set_ydata(FeasibleInput2)
		chosenpoint.set_xdata(U[0,i])
		chosenpoint.set_ydata(U[1,i])
		TimeText.set_text("t = " + str(t[i]))
		return feasible,cline,cline_dashed,chosenpoint,TimeText,


	# Init only required for blitting to give a clean slate.
	def init():
		ax1.plot([Bounds[0][0],Bounds[0][1]],[Bounds[1][0],Bounds[1][0]],'k--')
		ax1.plot([Bounds[0][0],Bounds[0][1]],[Bounds[1][1],Bounds[1][1]],'k--')
		ax1.plot([Bounds[0][0],Bounds[0][0]],[Bounds[1][0],Bounds[1][1]],'k--')
		ax1.plot([Bounds[0][1],Bounds[0][1]],[Bounds[1][0],Bounds[1][1]],'k--')
		Coefficient1,Coefficient2,Constraint1 = return_constraint_variables_muscle_activation_driven(t[0],X[:,0])
		AllowableBounds_x = np.array([U[0,0]-MaxStep,U[0,0]+MaxStep])
		AllowableBounds_y = np.array([U[1,0]-MaxStep,U[1,0]+MaxStep])

		if Coefficient1 == 0:
			LowerBound_x = max(Bounds[0][0],AllowbaleBounds_x[0])
			UpperBound_x = min(Bounds[0][1],AllowbaleBounds_x[1])
			LowerBound_y = Constraint1/Coefficient2
			UpperBound_y = Constraint1/Coefficient2
			FeasibleInput1 = (UpperBound_x-LowerBound_x)*np.random.rand(1000) + LowerBound_x
			FeasibleInput2 = np.array([Constraint1/Coefficient2]*1000)
		elif Coefficient2 == 0:
			LowerBound_x = Constraint1/Coefficient1
			UpperBound_x = Constraint1/Coefficient1
			LowerBound_y = max(Bounds[1][0],AllowableBounds_y[0])
			UpperBound_y = min(Bounds[1][1],AllowableBounds_y[1])
			FeasibleInput1 = np.array([Constraint1/Coefficient1]*1000)
			FeasibleInput2 = (UpperBound_y-LowerBound_y)*np.random.rand(1000) + LowerBound_y
		else:
			SortedAllowableBounds = np.sort([\
										(Constraint1-Coefficient2*AllowableBounds_y[0])/Coefficient1,\
										(Constraint1-Coefficient2*AllowableBounds_y[1])/Coefficient1\
										])
			SortedBounds = np.sort([(Constraint1-Coefficient2*Bounds[1][0])/Coefficient1,\
										(Constraint1-Coefficient2*Bounds[1][1])/Coefficient1])
			LowerBound_x = max(	Bounds[0][0],\
			 					SortedBounds[0],\
								AllowableBounds_x[0],\
								SortedAllowableBounds[0]\
							)
			UpperBound_x = min(	Bounds[0][1],\
			 					SortedBounds[1],\
								AllowableBounds_x[1],\
								SortedAllowableBounds[1]\
							)
			LowerBound_y,UpperBound_y = np.sort([\
										(Constraint1-Coefficient1*LowerBound_x)/Coefficient2,\
										(Constraint1-Coefficient1*UpperBound_x)/Coefficient2\
										])
			# if UpperBound < LowerBound: import ipdb; ipdb.set_trace()
			assert UpperBound_x >= LowerBound_x, "Error generating bounds. Not feasible!"
			FeasibleInput1 = (UpperBound_x-LowerBound_x)*np.random.rand(1000) + LowerBound_x
			FeasibleInput2 = np.array([Constraint1/Coefficient2 - (Coefficient1/Coefficient2)*el \
									for el in FeasibleInput1])
		# feasible = plt.Circle((U[:,0]),radius=MaxStep_Activation,Color='b',alpha=0.5)
		feasible = plt.Rectangle((LowerBound_x,LowerBound_y),\
									UpperBound_x-LowerBound_x,\
									UpperBound_y-LowerBound_y,\
									alpha = 0.5,
									Color = 'b')
		feasible.set_visible(False)
		if Coefficient2!=0:
			cline_dashed, = plt.plot(\
				ax1.get_xlim(),\
				np.array(list(map(lambda x: \
					(Constraint1-Coefficient1*x)/Coefficient2, ax1.get_xlim()))),\
				'k--',lw=1)
		else:
			cline_dashed, = plt.plot(\
				np.array(list(map(lambda y: \
					(Constraint1-Coefficient2*y)/Coefficient1, ax1.get_ylim()))),\
				ax1.get_ylim(),\
				'k--',lw=1)
		cline_dashed.set_visible(False)
		cline, = plt.plot(FeasibleInput1,FeasibleInput2,'b',lw=2)
		cline.set_visible(False)
		chosenpoint, = plt.plot(U[0,0],U[1,0],c='k',marker='o')
		chosenpoint.set_visible(False)
		TimeText = plt.text(0.75,0.75,"t = " + str(t[0]),fontsize=16)
		TimeText.set_visible(False)
		return feasible,cline,cline_dashed,chosenpoint,TimeText,

	ani = animation.FuncAnimation(fig, animate, np.arange(1, np.shape(X)[1],1), init_func=init,interval=1, blit=False)
	plt.show()
def plot_individual_constraint_versus_time_muscle_activation_driven(
		t,X,**kwargs):
	"""
	A⋅u₁ + B⋅u₂ = C

	Takes in Time (t - numpy.ndarray of shape (N,)), the state array (X - numpy.ndarray of shape (4,N)) and plots the constraint equation and its coefficients over time.

	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	**kwargs
	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	1) Return - must be a bool. Determines whether or not a function handle is returned. Default is False.

	"""
	import numpy as np
	import matplotlib.pyplot as plt

	assert np.shape(X)[0] == 8 and str(type(X)) == "<class 'numpy.ndarray'>", "X must be a (8,N) numpy.ndarray"

	Return = kwargs.get("Return",False)
	assert type(Return) == bool, "Return must be either True or False."

	DescriptiveTitle = "Plotting Coefficients/Constraints vs. Time\nMuscle Velocity Driven\n$A\cdot u_{1} + B\cdot u_{2} = C$"

	fig, (ax1,ax2,ax3) = plt.subplots(1,3,figsize=(12,5))
	plt.suptitle(DescriptiveTitle,Fontsize=20,y=0.95)
	plt.subplots_adjust(wspace = 0.4,top=0.8)

	"""
	A⋅u₁ + B⋅u₂ = C
	"""

	A,B,C = [],[],[]
	for i in range(np.shape(X)[1]):
		Coefficient1,Coefficient2,Constraint1 = return_constraint_variables_muscle_activation_driven(t[i],X[:,i])
		A.append(Coefficient1)
		B.append(Coefficient2)
		C.append(Constraint1)

	ax1.plot(t[:np.shape(X)[1]],A,'r',lw=2)
	ax1.spines['right'].set_visible(False)
	ax1.spines['top'].set_visible(False)
	ax1.set_ylabel(r"$1^{st}$ Coefficient")
	ax1.set_xlabel("Time (s)")

	ax2.plot(t[:np.shape(X)[1]],B,'b',lw=2)
	ax2.spines['right'].set_visible(False)
	ax2.spines['top'].set_visible(False)
	ax2.set_ylabel(r"$2^{nd}$ Coefficient")
	ax2.set_xticks(ax1.get_xticks())
	ax2.set_xticklabels([""]*len(ax1.get_xticks()))

	ax3.plot(t[:np.shape(X)[1]],C,'k',lw=2)
	ax3.spines['right'].set_visible(False)
	ax3.spines['top'].set_visible(False)
	ax3.set_ylabel("Constraint")
	ax3.set_xticks(ax1.get_xticks())
	ax3.set_xticklabels([""]*len(ax1.get_xticks()))

	if Return == True:
		return(fig)
	else:
		plt.show()
def plot_individual_coefficient2_versus_time_muscle_activation_driven(
	t,X,**kwargs):
	"""
	B = c2⋅c4⋅c10⋅R2(X)⋅KT_2(X)⋅FLV_2(X)

	Takes in Time (t - numpy.ndarray of shape (N,)), the state array (X - numpy.ndarray of shape (4,N)) and plots the 2nd Coefficient of the Constraint Equation over time as well as its components.

	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	**kwargs
	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	1) Return - must be a bool. Determines whether or not a function handle is returned. Default is False.

	"""
	import numpy as np
	import matplotlib.pyplot as plt

	assert np.shape(X)[0] == 8 and str(type(X)) == "<class 'numpy.ndarray'>", "X must be a (8,N) numpy.ndarray"

	Return = kwargs.get("Return",False)
	assert type(Return) == bool, "Return must be either True or False."

	fig, (ax1,ax2,ax3,ax4) = plt.subplots(1,4,figsize=(15,5))
	plt.subplots_adjust(top=0.8,bottom=0.1,left=0.1,right=0.975,wspace=0.4)
	DescriptiveTitle = "Plotting $2^{nd}$ Coefficient vs. Time\n$B=c_{2}c_{4}c_{10}R_{2}(\\vec{x}(t))K_{T,2}(\\vec{x}(t))F_{LV,2}(\\vec{x}(t))$"
	plt.suptitle(DescriptiveTitle,Fontsize=20,y=0.95)


	r2,kt_2,flv_2,B = [],[],[],[]
	for i in range(np.shape(X)[1]):
		_,Coefficient2,_ = return_constraint_variables_muscle_activation_driven(t[i],X[:,i])
		B.append(Coefficient2)
		r2.append(R2(X[:,i]))
		kt_2.append(KT_2(X[:,i]))
		flv_2.append(FLV_2(X[:,i]))

	ax1.plot(t[:np.shape(X)[1]],r2,'b',lw=2)
	ax1.spines['right'].set_visible(False)
	ax1.spines['top'].set_visible(False)
	ax1.set_ylabel(r"$R_{2}(\vec{x}(t))$")
	ax1.set_xlabel("Time (s)")

	ax2.plot(t[:np.shape(X)[1]],kt_2,'b',lw=2)
	ax2.spines['right'].set_visible(False)
	ax2.spines['top'].set_visible(False)
	ax2.set_ylabel(r"$K_{T,2}(\vec{x}(t))$")
	ax2.set_xticks(ax1.get_xticks())
	ax2.set_xticklabels([""]*len(ax1.get_xticks()))

	ax3.plot(t[:np.shape(X)[1]],flv_2,'b',lw=2)
	ax3.spines['right'].set_visible(False)
	ax3.spines['top'].set_visible(False)
	ax3.set_ylabel(r"$F_{LV,2}(\vec{x}(t))$")
	ax3.set_xticks(ax1.get_xticks())
	ax3.set_xticklabels([""]*len(ax1.get_xticks()))

	ax4.plot(t[:np.shape(X)[1]],B,'b',lw=2)
	ax4.spines['right'].set_visible(False)
	ax4.spines['top'].set_visible(False)
	ax4.set_ylabel(r"$2^{nd}$ Coefficient")
	ax4.set_xticks(ax1.get_xticks())
	ax4.set_xticklabels([""]*len(ax1.get_xticks()))

	if Return == True:
		return(fig)
	else:
		plt.show()
def plot_individual_coefficient1_versus_time_muscle_activation_driven(
	t,X,**kwargs):
	"""
	A = c2⋅c3⋅c6⋅R1(X)⋅KT_1(X)⋅FLV_1(X)

	Takes in Time (t - numpy.ndarray of shape (N,)), the state array (X - numpy.ndarray of shape (4,N)) and plots the 1st Coefficient of the Constraint Equation over time as well as its components.

	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	**kwargs
	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	1) Return - must be a bool. Determines whether or not a function handle is returned. Default is False.

	"""
	import numpy as np
	import matplotlib.pyplot as plt

	assert np.shape(X)[0] == 8 and str(type(X)) == "<class 'numpy.ndarray'>", "X must be a (8,N) numpy.ndarray"

	Return = kwargs.get("Return",False)
	assert type(Return) == bool, "Return must be either True or False."

	fig, (ax1,ax2,ax3,ax4) = plt.subplots(1,4,figsize=(15,5))
	plt.subplots_adjust(top=0.8,bottom=0.1,left=0.1,right=0.975,wspace=0.4)
	DescriptiveTitle = "Plotting $1^{st}$ Coefficient vs. Time\n$A=c_{2}c_{3}c_{6}R_{1}(\\vec{x}(t))K_{T,1}(\\vec{x}(t))F_{LV,1}(\\vec{x}(t))$"
	plt.suptitle(DescriptiveTitle,Fontsize=20,y=0.95)


	r1,kt_1,flv_1,A = [],[],[],[]
	for i in range(np.shape(X)[1]):
		Coefficient1,_,_ = return_constraint_variables_muscle_activation_driven(t[i],X[:,i])
		A.append(Coefficient1)
		r1.append(R1(X[:,i]))
		kt_1.append(KT_1(X[:,i]))
		flv_1.append(FLV_1(X[:,i]))

	ax1.plot(t[:np.shape(X)[1]],r1,'r',lw=2)
	ax1.spines['right'].set_visible(False)
	ax1.spines['top'].set_visible(False)
	ax1.set_ylabel(r"$R_{1}(\vec{x}(t))$")
	ax1.set_xlabel("Time (s)")

	ax2.plot(t[:np.shape(X)[1]],kt_1,'r',lw=2)
	ax2.spines['right'].set_visible(False)
	ax2.spines['top'].set_visible(False)
	ax2.set_ylabel(r"$K_{T,1}(\vec{x}(t))$")
	ax2.set_xticks(ax1.get_xticks())
	ax2.set_xticklabels([""]*len(ax1.get_xticks()))

	ax3.plot(t[:np.shape(X)[1]],flv_1,'r',lw=2)
	ax3.spines['right'].set_visible(False)
	ax3.spines['top'].set_visible(False)
	ax3.set_ylabel(r"$F_{LV,1}(\vec{x}(t))$")
	ax3.set_xticks(ax1.get_xticks())
	ax3.set_xticklabels([""]*len(ax1.get_xticks()))

	ax4.plot(t[:np.shape(X)[1]],A,'r',lw=2)
	ax4.spines['right'].set_visible(False)
	ax4.spines['top'].set_visible(False)
	ax4.set_ylabel(r"$1^{st}$ Coefficient")
	ax4.set_xticks(ax1.get_xticks())
	ax4.set_xticklabels([""]*len(ax1.get_xticks()))

	if Return == True:
		return(fig)
	else:
		plt.show()

def plot_states(t,X,**kwargs):
	"""
	Takes in a numpy.ndarray for time (t) of shape (N,) and the numpy.ndarray for the state space (X) of shape (M,N), where M is the number of states and N is the same length as time t. Returns a plot of the states.

	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	**kwargs
	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	1) Return - must be a bool. Determines if the function returns a function handle. Default is False.

	2) InputString - must be a string. Input to the DescriptiveTitle that can be used to personalize the title. Default is None.

	"""
	import numpy as np
	import matplotlib.pyplot as plt

	assert (np.shape(X)[0] in [2,4,8]) \
				and (np.shape(X)[1] == len(t)) \
					and (str(type(X)) == "<class 'numpy.ndarray'>"), \
			"X must be a (2,N), (4,N), or (8,N) numpy.ndarray, where N is the length of t."


	Return = kwargs.get("Return",False)
	assert type(Return)==bool, "Return must be either True or False."

	InputString = kwargs.get("InputString",None)
	assert InputString is None or type(InputString)==str, "InputString must either be None or a str."

	NumStates = np.shape(X)[0]
	X[:2,:] = 180*X[:2,:]/np.pi # converting to deg and deg/s
	NumRows = int(np.ceil(NumStates/5))
	if NumStates < 5:
		NumColumns = NumStates
	else:
		NumColumns = 5

	ColumnNumber = [el%5 for el in np.arange(0,NumStates,1)]
	RowNumber = [int(el/5) for el in np.arange(0,NumStates,1)]
	Units = ["(Deg)","(Deg/s)","(N)","(N)","(m)","(m)","(m/s)","(m/s)"]
	if InputString is None:
		DescriptiveTitle = "Plotting States vs. Time"
	else:
		assert type(InputString)==str, "InputString must be a string"
		DescriptiveTitle = InputString + " Driven"

	Figure = kwargs.get("Figure",None)
	assert (Figure is None) or \
				(	(type(Figure)==tuple) and \
					(str(type(Figure[0]))=="<class 'matplotlib.figure.Figure'>") and\
					(np.array([str(type(ax))=="<class 'matplotlib.axes._subplots.AxesSubplot'>" \
						for ax in Figure[1].flatten()]).all()) and \
					(Figure[1].shape == (NumRows,NumColumns))\
				),\
				 	"Figure can either be left blank (None) or it must be constructed from data that has the same shape as X."
	if Figure is None:
		fig, axes = plt.subplots(NumRows,NumColumns,figsize=(3*NumColumns,2*NumRows + 2))
		plt.subplots_adjust(top=0.85,bottom=0.15,left=0.075,right=0.975)
		plt.suptitle(DescriptiveTitle,Fontsize=20,y=0.975)
		if NumStates <=5:
			for j in range(NumStates):
				axes[ColumnNumber[j]].spines['right'].set_visible(False)
				axes[ColumnNumber[j]].spines['top'].set_visible(False)
				axes[ColumnNumber[j]].plot(t,X[j,:])
				if ColumnNumber[j]!=0:
					axes[ColumnNumber[j]].set_xticklabels(\
										[""]*len(axes[ColumnNumber[j]].get_xticks()))
				else:
					axes[ColumnNumber[j]].set_xlabel("Time (s)")
				axes[ColumnNumber[j]].set_title(r"$x_{" + str(j+1) + "}$ " + Units[j])
		else:
			for j in range(NumStates):
				axes[RowNumber[j],ColumnNumber[j]].spines['right'].set_visible(False)
				axes[RowNumber[j],ColumnNumber[j]].spines['top'].set_visible(False)
				axes[RowNumber[j],ColumnNumber[j]].plot(t,X[j,:])
				if not(RowNumber[j] == RowNumber[-1] and ColumnNumber[j]==0):
					axes[RowNumber[j],ColumnNumber[j]].set_xticklabels(\
										[""]*len(axes[RowNumber[j],ColumnNumber[j]].get_xticks()))
				else:
					axes[RowNumber[j],ColumnNumber[j]].set_xlabel("Time (s)")
				axes[RowNumber[j],ColumnNumber[j]].set_title(r"$x_{" + str(j+1) + "}$ "+ Units[j])
			if NumStates%5!=0:
				[fig.delaxes(axes[RowNumber[-1],el]) for el in range(ColumnNumber[-1]+1,5)]
	else:
		fig = Figure[0]
		axes = Figure[1]
		for i in range(NumStates):
			axes[RowNumber[i],ColumnNumber[i]].plot(t,X[i,:])
	X[:2,:] = np.pi*X[:2,:]/180
	if Return == True:
		return((fig,axes))
	else:
		plt.show()
def plot_inputs(t,U,**kwargs):
	"""
	Takes in a numpy.ndarray for time (t) of shape (N,) and the numpy.ndarray for the input (U) (NOT NECESSARILY THE SAME LENGTH AS t). Returns a plot of the states.

	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	**kwargs
	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	1) Return - must be a bool. Determines if the function returns a function handle. Default is False.

	2) InputString - must be a string. Input to the DescriptiveTitle that can be used to personalize the title. Default is None.

	"""
	import numpy as np
	import matplotlib.pyplot as plt

	assert (np.shape(U)[0] == 2) \
				and (np.shape(U)[1] == len(t)) \
					and (str(type(U)) == "<class 'numpy.ndarray'>"), \
			"X must be a (2,N) numpy.ndarray, where N is the length of t."

	Return = kwargs.get("Return",False)
	assert type(Return)==bool, "Return must be either True or False."

	InputString = kwargs.get("InputString",None)
	assert InputString is None or type(InputString)==str, "InputString must either be None or a str."

	if InputString is None:
		DescriptiveTitle = "Plotting Inputs vs. Time"
	else:
		assert type(InputString)==str, "InputString must be a string"
		DescriptiveTitle = InputString + " vs. Time"

	Figure = kwargs.get("Figure",None)
	assert (Figure is None) or \
				(	(type(Figure)==tuple) and \
					(str(type(Figure[0]))=="<class 'matplotlib.figure.Figure'>") and\
					(np.array([str(type(ax))=="<class 'matplotlib.axes._subplots.AxesSubplot'>" \
						for ax in Figure[1].flatten()]).all()) and \
					(Figure[1].shape == (2,))\
				),\
				 	"Figure can either be left blank (None) or it must be constructed from data that has the same shape as U."
	if Figure is None:
		fig, axes = plt.subplots(1,2,figsize=(13,5))
		plt.subplots_adjust(top=0.9,hspace=0.4,bottom=0.1,left=0.075,right=0.975)
		plt.suptitle(DescriptiveTitle,Fontsize=20,y=0.975)

		axes[0].plot(t,U[0,:],lw=2)
		axes[0].plot([-1,t[-1]+1],[0,0],'k--',lw=0.5)
		axes[0].set_xlim([t[0],t[-1]])
		axes[0].spines['right'].set_visible(False)
		axes[0].spines['top'].set_visible(False)
		axes[0].set_ylabel(r"$u_1$")
		axes[0].set_xlabel("Time (s)")

		axes[1].plot(t,U[1,:],lw=2)
		axes[1].plot([-1,t[-1]+1],[0,0],'k--',lw=0.5)
		axes[1].set_xlim([t[0],t[-1]])
		axes[1].spines['right'].set_visible(False)
		axes[1].spines['top'].set_visible(False)
		axes[1].set_ylabel(r"$u_2$")
		axes[1].set_xticks(axes[0].get_xticks())
		axes[1].set_xticklabels([""]*len(axes[0].get_xticks()))
	else:
		fig = Figure[0]
		axes = Figure[1]
		axes[0].plot(t,U[0,:],lw=2)
		axes[1].plot(t,U[1,:],lw=2)

	if Return == True:
		return((fig,axes))
	else:
		plt.show()
def plot_l_m_comparison(t,X,**kwargs):

	"""
	Takes in a numpy.ndarray for time (t) of shape (N,) and the numpy.ndarray for the input (U) (NOT NECESSARILY THE SAME LENGTH AS t). Returns a plot of the states.

	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	**kwargs
	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	1) Return - must be a bool. Determines if the function returns a function handle. Default is False.

	2) InputString - must be a string. Input to the DescriptiveTitle that can be used to personalize the title. Default is None.

	"""
	import numpy as np
	import matplotlib.pyplot as plt

	assert (np.shape(X)[0] >= 2) \
				and (np.shape(X)[1] == len(t)) \
					and (str(type(X)) == "<class 'numpy.ndarray'>"), \
			"X must be a (M,N) numpy.ndarray, where M is greater than or equal to 2 and N is the length of t."

	Return = kwargs.get("Return",False)
	assert type(Return)==bool, "Return must be either True or False."

	InputString = kwargs.get("InputString",None)
	assert InputString is None or type(InputString)==str, "InputString must either be None or a str."
	if InputString is None:
		DescriptiveTitle = "Muscle vs. Musculotendon Lengths"
	else:
		DescriptiveTitle = "Muscle vs. Musculotendon Lengths\n" + InputString + " Driven"

	L_m = kwargs.get("MuscleLengths",None)
	assert (L_m is None) or (str(type(L_m))=="<class 'numpy.ndarray'>" and np.shape(L_m)==(2,len(t))), "L_m must either be a numpy.ndarray of size (2,N) or left as None (Default)."

	V_m = kwargs.get("MuscleVelocities",None)
	assert (V_m is None) or (str(type(V_m))=="<class 'numpy.ndarray'>" and np.shape(V_m)==(2,len(t))), "V_m must either be a numpy.ndarray of size (2,N) or left as None (Default)."

	assert L_m is not None or V_m is not None, "Error! Need to input some length/velocity measurement for the muscles."

	if L_m is None:
		"""
		This is for the muscle velocity driven controller. These values of initial muscle length are estimates taken to be the optimal muscle lengths. We will need to run some sensitivity analysis to ensure that this does not drastically effect the deviations from the MTU estimate.
		"""
		l_m1 = integrate.cumtrapz(V_m[0,:],t,initial = 0) + np.ones(len(t))*lo1
		l_m2 = integrate.cumtrapz(V_m[1,:],t,initial = 0) + np.ones(len(t))*lo2
	else:
		l_m1 = L_m[0,:]
		l_m2 = L_m[1,:]

	"""
	Note: X must be transposed in order to run through map()
	"""
	Figure = kwargs.get("Figure",None)
	assert (Figure is None) or \
				(	(type(Figure)==tuple) and \
					(str(type(Figure[0]))=="<class 'matplotlib.figure.Figure'>") and\
					(np.array([str(type(ax))=="<class 'matplotlib.axes._subplots.AxesSubplot'>" \
						for ax in Figure[1].flatten()]).all()) and \
					(Figure[1].shape == (2,2))\
				),\
				 	"Figure can either be left blank (None) or it must be constructed from data that has the same shape as X."

	if Figure is None:
		fig, axes = plt.subplots(2,2,figsize = (14,7))
		plt.suptitle(DescriptiveTitle,Fontsize=20,y=0.975)

		axes[0,0].plot(t,integrate.cumtrapz(np.array(list(map(lambda X: v_MTU1(X),X.T))),\
							t,initial=0) + np.ones(len(t))*l_m1[0], '0.70')
		axes[0,0].plot(t,l_m1)
		axes[0,0].set_ylabel(r"$l_{m,1}/l_{MTU,1}$ (m)")
		axes[0,0].set_xlabel("Time (s)")

		axes[0,1].plot(t,l_m1-integrate.cumtrapz(np.array(list(map(lambda X: v_MTU1(X),X.T))),\
							t,initial=0) - np.ones(len(t))*l_m1[0])
		axes[0,1].set_ylabel("Error (m)")
		axes[0,1].set_xlabel("Time (s)")

		axes[1,0].plot(t,integrate.cumtrapz(np.array(list(map(lambda X: v_MTU2(X),X.T))),\
							t,initial=0) + np.ones(len(t))*l_m2[0], '0.70')
		axes[1,0].plot(t,l_m2)
		axes[1,0].set_ylabel(r"$l_{m,2}/l_{MTU,2}$ (m)")
		axes[1,0].set_xlabel("Time (s)")

		axes[1,1].plot(t,l_m2-integrate.cumtrapz(np.array(list(map(lambda X: v_MTU2(X),X.T))),\
							t,initial=0) - np.ones(len(t))*l_m2[0])
		axes[1,1].set_ylabel("Error (m)")
		axes[1,1].set_xlabel("Time (s)")
	else:
		fig = Figure[0]
		axes = Figure[1]
		axes[0,0].plot(t,integrate.cumtrapz(np.array(list(map(lambda X: v_MTU1(X),X.T))),\
							t,initial=0) + np.ones(len(t))*l_m1[0], '0.70')
		axes[0,0].plot(t,l_m1)

		axes[0,1].plot(t,l_m1-integrate.cumtrapz(np.array(list(map(lambda X: v_MTU1(X),X.T))),\
							t,initial=0) - np.ones(len(t))*l_m1[0])

		axes[1,0].plot(t,integrate.cumtrapz(np.array(list(map(lambda X: v_MTU2(X),X.T))),\
							t,initial=0) + np.ones(len(t))*l_m2[0], '0.70')
		axes[1,0].plot(t,l_m2)

		axes[1,1].plot(t,l_m2-integrate.cumtrapz(np.array(list(map(lambda X: v_MTU2(X),X.T))),\
							t,initial=0) - np.ones(len(t))*l_m2[0])

	if Return == True:
		return((fig,axes))
	else:
		plt.show()
def plot_N_sim_MAT(t,TotalX,TotalU,**kwargs):
	Return = kwargs.get("Return",False)
	assert type(Return) == bool, "Return should either be True or False"

	fig1 = plt.figure(figsize = (9,7))
	fig1_title = "Underdetermined Forced-Pendulum Example"
	plt.title(fig1_title,fontsize=16,color='gray')
	statusbar = dsb(0,np.shape(TotalX)[0],title=(plot_N_sim_MAT.__name__ + " (" + fig1_title +")"))
	for j in range(np.shape(TotalX)[0]):
		plt.plot(t,(TotalX[j,0,:])*180/np.pi,'0.70',lw=2)
		statusbar.update(j)
	plt.plot(np.linspace(0,t[-1],1001),\
			(r(np.linspace(0,t[-1],1001)))*180/np.pi,\
				'r')
	plt.xlabel("Time (s)")
	plt.ylabel("Desired Measure (Deg)")

	fig2 = plt.figure(figsize = (9,7))
	fig2_title = "Error vs. Time"
	plt.title(fig2_title)
	statusbar.reset(title=(plot_N_sim_MAT.__name__ + " (" + fig2_title +")"))
	for j in range(np.shape(TotalX)[0]):
		plt.plot(t, (r(t)-TotalX[j,0,:])*180/np.pi,color='0.70')
		statusbar.update(j)
	plt.xlabel("Time (s)")
	plt.ylabel("Error (Deg)")

	statusbar.reset(
		title=(
			plot_N_sim_MAT.__name__
			+ " (Plotting States, Inputs, and Muscle Length Comparisons)"
			)
		)
	for j in range(np.shape(TotalX)[0]):
		if j == 0:
			fig3 = plot_states(t,TotalX[j],Return=True,InputString = "Muscle Activations")
			fig4 = plot_inputs(t,TotalU[j],Return=True,InputString = "Muscle Activations")
			fig5 = plot_l_m_comparison(t,TotalX[j],MuscleLengths = TotalX[j,4:6,:],Return=True,InputString = "Muscle Activation")
		else:
			fig3 = plot_states(t,TotalX[j],Return=True,InputString = "Muscle Activations",\
									Figure=fig3)
			fig4 = plot_inputs(t,TotalU[j],Return=True,InputString = "Muscle Activations", \
									Figure = fig4)
			fig5 = plot_l_m_comparison(t,TotalX[j],MuscleLengths = TotalX[j,4:6,:],Return=True,\
											InputString = "Muscle Activation", Figure = fig5)
		statusbar.update(j)
	if Return == True:
		return([fig1,fig2,fig3,fig4,fig5])
	else:
		plt.show()

def save_figures(BaseFileName,figs):
	import os.path
	from matplotlib.backends.backend_pdf import PdfPages
	i = 1
	FileName = BaseFileName + "_" + "{:0>2d}".format(i) + ".pdf"
	if os.path.exists(FileName) == True:
		while os.path.exists(FileName) == True:
			i += 1
			FileName = BaseFileName + "_" + "{:0>2d}".format(i) + ".pdf"
	PDFFile = PdfPages(FileName)
	if len(figs)==1:
		PDFFile.savefig(figs[0])
	else:
		[PDFFile.savefig(fig) for fig in figs]
	PDFFile.close()

N = N_seconds*10000 + 1
t = np.linspace(0,N_seconds,N)
dt = t[1]-t[0]

TotalX,TotalU = run_N_sim_MAT(NumberOfTrials=3)

# plot_N_sim_MAT(t,TotalX,TotalU,Return=False)
figs = plot_N_sim_MAT(t,TotalX,TotalU,Return=True)
figs=[manager.canvas.figure
         for manager in matplotlib._pylab_helpers.Gcf.get_all_fig_managers()]
save_figures("1DOF_2DOA_Minimum_Activation_Transition",figs)
plt.close('all')
