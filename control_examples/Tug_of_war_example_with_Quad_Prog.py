import numpy as np
import matplotlib.pyplot as plt
from danpy.sb import dsb
from scipy.optimize import fsolve
import cvxpy as cp

N_seconds = 20
dt = 0.01
t = np.arange(0,N_seconds + dt, dt)

us1 = 1
us2 = 1

K = 1
B = 1
M = 1

k1 = 10
k2 = 10

k_0 = 1
k_1 = 1

def k1_func(x1,u1):
    # return(k1)
    assert (x1 + (u1-us1))>=-1e-6
    return(k1)

def k2_func(x1,u2):
    # return(k2)
    assert (x1 - (u2-us2))>=-1e-6
    return(k2)


def dx1_dt(x1,x2,u1,u2):
    return(x2)

def dx2_dt(x1,x2,u1,u2):
    return(
        -K/M*x1
        - B/M*x2
        - k1*(x1 + (u1-us1))/M
        - k2*(x1 - (u2-us2))/M
    )

X = np.zeros((2,len(t)))
X[0,0] = -0.1
X[1,0] = 0
U = np.zeros((2,len(t)))
statusbar = dsb(0,len(t)-1,title="What are you doing?...")
P = np.array([[1,0],[0,1]])
q = np.array([[0],[0]])
G = -np.eye(2)
A = np.array([k1,-k2])
for i in range(len(t)-1):
    h = -np.array([us1,us2]) + np.array([1,-1])*X[0,i]
    u = cp.Variable(2)
    b = np.array(
        -(K+k1+k2)*X[0,i]
        - B*X[1,i]
        + k1*us1
        - k2*us2
        + M*k_0*X[0,i]
        + M*k_1*X[1,i]
    )
    prob = cp.Problem(cp.Minimize((1/2)*cp.quad_form(u, P) + q.T@u),
                 [G@u <= h,
                  A@u == b])
    prob.solve()
    u1,u2 = u.value
    # u1p,u2p = (1/np.sqrt(k1**2+k2**2))*np.array([[k1],[-k2]])*(
    #     -(K+k1+k2)*X[0,i]
    #     - B*X[1,i]
    #     + k1*us1
    #     - k2*us2
    #     + M*k_0*X[0,i]
    #     + M*k_1*X[1,i]
    #     )
    # c1 = fsolve(lambda c:u1p+c*k2-us1+X[0,i],0.1)
    # c2 = fsolve(lambda c:u2p+c*k1-us2-X[0,i],0.1)
    # if u1p+c1*k2-us1+X[0,i]>=-1e-6 and u2p+c1*k1-us2-X[0,i]>=-1e-6:
    #     u1 = u1p+c1*k2
    #     u2 = u2p+c1*k1
    # elif u1p+c2*k2-us1+X[0,i]>=-1e-6 and u2p+c2*k1-us2-X[0,i]>=-1e-6:
    #     u1 = u1p+c2*k2
    #     u2 = u2p+c2*k1
    # else:
    #     import ipdb; ipdb.set_trace()
    # if X[0,i]+u1<us1:
    #     u1 = us1-X[0,i]
    #     u2 = (-1/k2)*(
    #         -k1*u1
    #         - (K+k1+k2)*X[0,i]
    #         - B*X[1,i]
    #         + k1*us1
    #         - k2*us2
    #         + M*k_0*X[0,i]
    #         + M*k_1*X[1,i]
    #     )
    # elif u2-X[0,i]<us2:
    #     u2 = us2+X[0,i]
    #     u1 = (1/k1)*(
    #         -k2*u2
    #         - (K+k1+k2)*X[0,i]
    #         - B*X[1,i]
    #         + k1*us1
    #         - k2*us2
    #         + M*k_0*X[0,i]
    #         + M*k_1*X[1,i]
    #     )
    # while u1+X[0,i]<us1 or u2-X[0,i]<us2:
    #     u1=u1+0.1
    #     u2=u2+k1/k2*0.1
    U[:,i] = [u1,u2]
    # u1,u2 = 0,0
    X[0,i+1]=X[0,i]+dx1_dt(X[0,i],X[1,i],u1,u2)*dt
    X[1,i+1]=X[1,i]+dx2_dt(X[0,i],X[1,i],u1,u2)*dt
    statusbar.update(i)

plt.plot(t,X[0,:],'r')
plt.plot(t,X[1,:],'b')
plt.show()
