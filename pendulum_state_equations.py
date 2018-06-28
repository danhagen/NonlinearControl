from muscle_params_BIC_TRI_5_6 import *
import numpy as np

def v_MTU1(X):
	return(np.sign(-R1(X))*X[1]*np.sqrt(dR1_dx1(X)**2 + R1(X)**2)) # NOT NORMALIZED (in m/s)
def a_MTU1(X):
	return(np.sign(-R1(X))*(dX2_dt(X)*np.sqrt(dR1_dx1(X)**2 + R1(X)**2) \
				+ (X[1]**2)*dR1_dx1(X)*(d2R1_dx12(X) + R1(X))/np.sqrt(dR1_dx1(X)**2 + R1(X)**2)))

def v_MTU2(X):
	return(np.sign(-R2(X))*X[1]*np.sqrt(dR2_dx1(X)**2 + R2(X)**2)) # NOT NORMALIZED (in m/s)
def a_MTU2(X):
	return(np.sign(-R2(X))*(dX2_dt(X)*np.sqrt(dR2_dx1(X)**2 + R2(X)**2) \
				+ (X[1]**2)*dR2_dx1(X)*(d2R2_dx12(X) + R2(X))/np.sqrt(dR2_dx1(X)**2 + R2(X)**2)))


# g,L = 9.80, 0.45 #m/s², m
g,L = 0, 0.45 #m/s², m REMOVING GRAVITY
M = 1.6 # kg

c1 = -(3*g)/(2*L)
c2 = 3/(M*L**2)
c3 = np.cos(α1)
c4 = np.cos(α2)
c5 = np.cos(α1)/m1
c6 = F_MAX1*np.cos(α1)**2/m1
c7 = F_MAX1*bm1*np.cos(α1)**2/(m1*lo1)
c8 = np.tan(α1)**2
c9 = np.cos(α2)/m2
c10 = F_MAX2*np.cos(α2)**2/m2
c11 = F_MAX2*bm2*np.cos(α2)**2/(m2*lo2)
c12 = np.tan(α2)**2

def dX1_dt(X):
	return(X[1])
def d2X1_dt2(X):
	return(dX2_dt(X))
def dX2_dt(X,U=None):
	if U is None:
		return(c1*np.sin(X[0]) + c2*R1(X)*X[2] + c2*R2(X)*X[3])
	else:
		return(c1*np.sin(X[0]) + c2*R1(X)*U[0] + c2*R2(X)*U[1])
def d2X2_dt2(X):
	return(c1*np.cos(X[0])*dX1_dt(X) + c2*dR1_dx1(X)*dX1_dt(X)*X[2] + c2*R1(X)*dX3_dt(X)\
			+ c2*dR2_dx1(X)*dX1_dt(X)*X[3] + c2*R2(X)*dX4_dt(X))
def dX3_dt(X,U=None):
	if U is None:
		return(KT_1(X)*(v_MTU1(X) - c3*X[6]))
	else:
		return(KT_1(X)*(v_MTU1(X) - c3*U[0]))
def dX4_dt(X,U=None):
	if U is None:
		return(KT_2(X)*(v_MTU2(X) - c4*X[7]))
	else:
		return(KT_2(X)*(v_MTU2(X) - c4*U[1]))
def dX5_dt(X):
	return(X[6])
def dX6_dt(X):
	return(X[7])
def dX7_dt(X,U):
	return(c5*X[2] - c6*F_PE1_1(X) - c7*X[6] + c8*X[6]**2/X[4] - c6*FLV_1(X)*U[0])
def dX8_dt(X,U):
	return(c9*X[3] - c10*F_PE1_2(X) - c11*X[7] + c12*X[7]**2/X[5] - c10*FLV_2(X)*U[1])
