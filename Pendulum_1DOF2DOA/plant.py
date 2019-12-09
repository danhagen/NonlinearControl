import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from danpy.sb import dsb
from danpy.useful_functions import save_figures,is_number
from scipy import signal
import numdifftools as nd
import scipy as sp
from params import *

def LP_filt(filter_length, x):
	"""
	Finite Impulse Response (FIR) Moving Average (MA) Low-Pass Filter
	"""
	b=np.ones(filter_length,)/(filter_length) #Finite Impulse Response (FIR) Moving Average (MA) filter with one second filter length
	a=1
	y = signal.filtfilt(b, a, x)
	return y

class plant_pendulum_1DOF2DOF:
    def __init__(self,**params):
        self.Ij = params.get("Joint Inertia", 1.15e-2) # kg⋅m²
        is_number(self.Ij,"Joint Inertia",default=1.15e-2)

        self.bj = params.get("Joint Damping", 0.001) # N⋅s⋅m⁻¹
        is_number(self.bj,"Joint Damping",default=0.001)

        self.mj = params.get("Joint Mass", 0.541) # kg
        is_number(self.mj,"Joint Mass",default=0.541)

        self.rj = params.get("Joint Moment Arm", 0.05) # m
        is_number(self.rj,"Joint Moment Arm",default=0.05)

        self.Lcm = params.get("Link Center of Mass", 0.085) # m
        is_number(self.Lcm,"Link Center of Mass",default=0.085)

        self.L = params.get("Link Length", 0.3) # m
        is_number(self.L,"Link Length",default=0.3)

        self.Jm = params.get("Motor Inertia", 6.6e-5) # kg⋅m²
        is_number(self.Jm,"Motor Inertia",default=6.6e-5)

        self.bm = params.get("Motor Damping", 0.00462) # N⋅s⋅m⁻¹
        is_number(self.bm,"Motor Damping",default=0.00462)

        self.rm = params.get("Motor Moment Arm", 0.01) # m
        is_number(self.rm,"Motor Moment Arm",default=0.01)

        self.k_spr = params.get("Spring Stiffness Coefficient",1) # N
        is_number(self.k_spr,"",default=1)

        self.b_spr = params.get("Spring Shape Coefficient",100) # unit-less
        is_number(self.b_spr,"",default=1)

        self.simulationDuration = params.get("Simulation Duration", 1000)
        is_number(self.simulationDuration,"Simulation Duration")

        self.dt = params.get("dt", 0.01)
        is_number(self.dt,"dt")

        self.offset = params.get("offset", 0)
        is_number(self.offset,"offset")

        self.Lf4h0_list = []
        self.Lf2hs_list = []

        self.df2dx1_list = []
        self.df2dx2_list = []
        self.df2dx3_list = []
        self.df2dx5_list = []
        self.vs_list = []
    def C(self,X):
        """
        Returns zero until the effects are quantified
        """
        return(
            0
        )
    def dCdx1(self,X):
        return(0)
    def d2Cdx12(self,X):
        return(0)
    def d2Cdx1x2(self,X):
        return(0)
    def dCdx2(self,X):
        return(0)
    def d2Cdx22(self,X):
        return(0)

    def update_state_variables(self,X):

        #>>>> State functions

        self.f1 = self.f1_func(X)
        self.f2 = self.f2_func(X)
        self.f3 = self.f3_func(X)
        self.f4 = self.f4_func(X)
        self.f5 = self.f5_func(X)
        self.f6 = self.f6_func(X)

        #>>>> State functions first gradient

        # self.df1dx1 = 0
        self.df1dx2 = 1
        # self.df1dx3 = 0
        # self.df1dx4 = 0
        # self.df1dx5 = 0
        # self.df1dx6 = 0

        self.df2dx1 = self.df2dx1_func(X)
        self.df2dx1_list.append(self.df2dx1)
        self.df2dx2 = self.df2dx2_func(X)
        self.df2dx2_list.append(self.df2dx2)
        self.df2dx3 = self.df2dx3_func(X)
        self.df2dx3_list.append(self.df2dx3)
        # self.df2dx4 = 0
        self.df2dx5 = self.df2dx5_func(X)
        self.df2dx5_list.append(self.df2dx5)
        # self.df2dx6 = 0

        # self.df3dx1 = 0
        # self.df3dx2 = 0
        # self.df3dx3 = 0
        self.df3dx4 = 1
        # self.df3dx5 = 0
        # self.df3dx6 = 0

        # self.df4dx1 = N/A
        # self.df4dx2 = N/A
        # self.df4dx3 = N/A
        # self.df4dx4 = N/A
        # self.df4dx5 = N/A
        # self.df4dx6 = N/A

        # self.df5dx1 = 0
        # self.df5dx2 = 0
        # self.df5dx3 = 0
        # self.df5dx4 = 0
        # self.df5dx5 = 0
        self.df5dx6 = 1

        # self.df6dx1 = N/A
        # self.df6dx2 = N/A
        # self.df6dx3 = N/A
        # self.df6dx4 = N/A
        # self.df6dx5 = N/A
        # self.df6dx6 = N/A

        #>>>> State functions second gradient

        self.d2f2dx12 = self.d2f2dx12_func(X)
        self.d2f2dx1x2 = self.d2f2dx1x2_func(X)
        self.d2f2dx1x3 = self.d2f2dx1x3_func(X)
        self.d2f2dx1x5 = self.d2f2dx1x5_func(X)

        self.d2f2dx22 = self.d2f2dx22_func(X)

        self.d2f2dx32 = self.d2f2dx32_func(X)

        self.d2f2dx52 = self.d2f2dx52_func(X)

    # def motor_coupling_function(self,X,motorNumber):
    #     return(
    #         self.rm*self.k_spr*(
    #             np.exp(
    #                 self.b_spr*(
    #                     self.rm*X[2+2*(motorNumber-1)]
    #                     + ((1.5-motorNumber)/0.5)*self.rj*X[0]
    #                 )
    #             )
    #             -1
    #         )
    #     )
    def tendon_1_FL_func(self,X):
        return(
            self.k_spr*(
                np.exp(self.b_spr*(self.rm*X[2]+self.rj*X[0]))
                - 1
            )
        )
    def tendon_2_FL_func(self,X):
        return(
            self.k_spr*(
                np.exp(self.b_spr*(self.rm*X[2]-self.rj*X[0]))
                - 1
            )
        )

    def f1_func(self,X):
        return(X[1])

    def f2_func(self,X):
        return(
            (
                -self.C(X) # Coriolis and centrifugal torques (zero)
                - self.bj*X[1] # damping torque
                - self.Lcm*self.mj*gr*np.cos(X[0]) # gravitational torque
                - self.rj*self.k_spr * (
                    np.exp(self.b_spr*(self.rm*X[2] + self.rj*X[0]))
                    - np.exp(self.b_spr*(self.rm*X[4] - self.rj*X[0]))
                ) # total coupling torque between motors and joint
            )/self.Ij
        )
    def df2dx1_func(self,X):
        return(
            (
                -self.dCdx1(X) # Coriolis and centrifugal torques (zero)
                + self.Lcm*self.mj*gr*np.sin(X[0]) # gravitational torque
                - (self.rj**2)*self.k_spr*self.b_spr * (
                    np.exp(self.b_spr*(self.rm*X[2] + self.rj*X[0]))
                    + np.exp(self.b_spr*(self.rm*X[4] - self.rj*X[0]))
                ) # total coupling torque between motors and joint
            )/self.Ij
        )
    def d2f2dx12_func(self,X):
        return(
            (
                -self.d2Cdx12(X) # Coriolis and centrifugal torques (zero)
                + self.Lcm*self.mj*gr*np.cos(X[0]) # gravitational torque
                - (self.rj**3)*self.k_spr*(self.b_spr**2) * (
                    np.exp(self.b_spr*(self.rm*X[2] + self.rj*X[0]))
                    - np.exp(self.b_spr*(self.rm*X[4] - self.rj*X[0]))
                ) # total coupling torque between motors and joint
            )/self.Ij
        )
    def d2f2dx1x2_func(self,X):
        return(
            (
                -self.d2Cdx1x2(X) # Coriolis and centrifugal torques (zero)
            )/self.Ij
        )
    def d2f2dx1x3_func(self,X):
        """
        This is equivalently -dSda/Ij
        """
        return(
            -(self.rj**2)*self.rm*self.k_spr*(self.b_spr**2) * (
                np.exp(self.b_spr*(self.rm*X[2] + self.rj*X[0]))
            ) / self.Ij
        )
    def d2f2dx1x5_func(self,X):
        """
        This is equivalently dSdb/Ij
        """
        return(
            -(self.rj**2)*self.rm*self.k_spr*(self.b_spr**2) * (
                np.exp(self.b_spr*(self.rm*X[4] - self.rj*X[0]))
            ) / self.Ij
        )
    def df2dx2_func(self,X):
        return(
            (
                -self.dCdx2(X) # Coriolis and centrifugal torques (zero)
                - self.bj # damping torque
            )/self.Ij
        )
    def d2f2dx22_func(self,X):
        return(
            (
                -self.d2Cdx22(X) # Coriolis and centrifugal torques (zero)
            )/self.Ij
        )
    def df2dx3_func(self,X):
        """
        Equivalently, this is the negative value of -Q_{11}/Ij
        """
        return(
            -self.rj*self.rm*self.k_spr*self.b_spr * (
                np.exp(self.b_spr*(self.rm*X[2] + self.rj*X[0]))
            ) / self.Ij
        )
    def d2f2dx32_func(self,X):
        return(
            -self.rj*(self.rm**2)*self.k_spr*(self.b_spr**2) * (
                np.exp(self.b_spr*(self.rm*X[2] + self.rj*X[0]))
            ) / self.Ij
        )
    def df2dx5_func(self,X):
        """
        Equivalently, this is the negative value of -Q_{12}/Ij
        """
        return(
            self.rj*self.rm*self.k_spr*self.b_spr * (
                np.exp(self.b_spr*(self.rm*X[4] - self.rj*X[0]))
            ) / self.Ij
        )
    def d2f2dx52_func(self,X):
        return(
            self.rj*(self.rm**2)*self.k_spr*(self.b_spr**2) * (
                np.exp(self.b_spr*(self.rm*X[4] - self.rj*X[0]))
            ) / self.Ij
        )

    def f3_func(self,X):
        return(X[3])

    def f4_func(self,X):
        return(
            (
                -self.bm*X[3]
                - self.rm*self.k_spr*(
                    np.exp(self.b_spr*(self.rm*X[2] + self.rj*X[0]))
                    -1
                )
            )/self.Jm
        )

    def f5_func(self,X):
        return(X[5])

    def f6_func(self,X):
        return(
            (
                -self.bm*X[5]
                - self.rm*self.k_spr*(
                    np.exp(self.b_spr*(self.rm*X[4] - self.rj*X[0]))
                    -1
                )
            )/self.Jm
        )

    def f(self,X):
        result = np.zeros((6,1))
        result[0,0] = self.f1_func(X)
        result[1,0] = self.f2_func(X)
        result[2,0] = self.f3_func(X)
        result[3,0] = self.f4_func(X)
        result[4,0] = self.f5_func(X)
        result[5,0] = self.f6_func(X)
        return(result)
    def g(self,X):
        result = np.matrix(np.zeros((6,2)))
        result[3,0] = 1/self.Jm
        result[5,1] = 1/self.Jm
        return(result)
    def h(self,X):
        result = np.zeros((2,))
        result[0] = X[0]
        result[1] = (self.rj**2)*self.k_spr*self.b_spr*(
            np.exp(self.b_spr*(self.rm*X[2] + self.rj*X[0]))
            + np.exp(self.b_spr*(self.rm*X[4] - self.rj*X[0]))
        )
        return(result)
    def create_derivative_functions(self):
        grad_f2_func = nd.Gradient(self.f2_func)
        hess_f2_func = nd.Hessian(self.f2_func)
        jac_f_func = nd.Jacobian(self.f)
        out = {
            "Gradient of f2" : grad_f2_func,
            "Hessian of f2" : hess_f2_func,
            "Jacobian of f" : jac_f_func
        }
        return(out)

    def forward_simulation(self,Time,X_o,U=None):
        assert len(X_o)==6, "X_o must have 6 elements, not " + str(len(X_o)) + "."
        dt = Time[1]-Time[0]
        if U is None:
            U = np.zeros((2,len(Time)-1))
        else:
            assert np.shape(U)==(2,len(Time)-1), "U must be either None (default) of have shape (2,len(Time)-1), not " + str(np.shape(U)) + "."
        X = np.zeros((6,len(Time)))
        Y = np.zeros((2,len(Time)))
        X[:,0] = X_o
        Y[:,0] = self.h(X[:,0])
        # self.update_state_variables(X_o)
        # deriv_funcs = self.create_derivative_functions()
        # self.grad_f2_func = deriv_funcs["Gradient of f2"]
        # self.hess_f2_func = deriv_funcs["Hessian of f2"]
        # self.jac_f_func = deriv_funcs["Jacobian of f"]
        statusbar=dsb(0,len(Time)-1,title="Forward Simulation (Custom)")
        for i in range(len(Time)-1):
            X[:,i+1] = (
                X[:,i]
                + dt*(
                    self.f(X[:,i])
                    + self.g(X[:,i])@U[:,np.newaxis,i]
                ).T
            )
            Y[:,i+1] = self.h(X[:,i+1])
            # self.update_state_variables(X[:,i+1])
            statusbar.update(i)
        return(X,U,Y)

    def h0(self,X):
        return(X[0])
    def Lfh0(self,X):
        return(X[1])
    def Lf2h0(self,X):
        return(self.f2_func(X))
    # def Lf3h0(self,X):
    #     return(
    #         (self.grad_f2_func(X)@self.f(X))[0]
    #     )
    def Lf3h0(self,X):
        return(
            self.df2dx1_func(X)*self.f1_func(X)
            + self.df2dx2_func(X)*self.f2_func(X)
            + self.df2dx3_func(X)*self.f3_func(X)
            + self.df2dx5_func(X)*self.f5_func(X)
        )
    # def Lf4h0(self,X):
    #     return(
    #         (
    #             (
    #                 self.grad_f2_func(X)@self.jac_f_func(X)
    #                 + (self.hess_f2_func(X)@self.f(X)).T
    #             )@self.f(X)
    #         )[0,0]
    #     )
    def Lf4h0(self,X):
        return(
            (
                self.d2f2dx12_func(X)*self.f1_func(X)
                + self.d2f2dx1x2_func(X)*self.f2_func(X)
                + self.df2dx2_func(X)*self.df2dx1_func(X)
                + self.d2f2dx1x3_func(X)*self.f3_func(X)
                + self.d2f2dx1x5_func(X)*self.f5_func(X)
            ) * self.f1_func(X)
            + (
                self.d2f2dx1x2_func(X)*self.f1_func(X)
                + self.df2dx1_func(X)
                + self.d2f2dx22_func(X)*self.f2_func(X)
                + (self.df2dx2_func(X)**2)
            ) * self.f2_func(X)
            + (
                self.d2f2dx1x3_func(X)*self.f1_func(X)
                + self.df2dx2_func(X)*self.df2dx3_func(X)
                + self.d2f2dx32_func(X)*self.f3_func(X)
            ) * self.f3_func(X)
            + (
                self.df2dx3_func(X)
            ) * self.f4_func(X)
            + (
                self.d2f2dx1x5_func(X)*self.f1_func(X)
                + self.df2dx2_func(X)*self.df2dx5_func(X)
                + self.d2f2dx52_func(X)*self.f5_func(X)
            ) * self.f5_func(X)
            + (
                self.df2dx5_func(X)
            ) * self.f6_func(X)
        )

    def hs(self,X):
        return(
            (self.rj**2)*self.k_spr*self.b_spr*(
                np.exp(self.b_spr*(self.rm*X[2] + self.rj*X[0]))
                + np.exp(self.b_spr*(self.rm*X[4] - self.rj*X[0]))
            )
        )
    def Lfhs(self,X):
        return(
            (self.rj**2)*self.k_spr*(self.b_spr**2)*(
                (self.rj*self.f1_func(X) + self.rm*self.f3_func(X))*(
                    np.exp(self.b_spr*(self.rm*X[2] + self.rj*X[0]))
                )
                - (self.rj*self.f1_func(X) - self.rm*self.f5_func(X))*(
                    np.exp(self.b_spr*(self.rm*X[4] - self.rj*X[0]))
                )
            )
        )
    def Lf2hs(self,X):
        return(
            (self.rj**2)*self.k_spr*(self.b_spr**2)*(
                (
                    self.b_spr*(self.rj*self.f1_func(X) + self.rm*self.f3_func(X))**2
                    + (self.rj*self.f2_func(X) + self.rm*self.f4_func(X))
                ) * np.exp(self.b_spr*(self.rm*X[2] + self.rj*X[0]))

                + (
                    self.b_spr*(self.rj*self.f1_func(X) - self.rm*self.f5_func(X))**2
                    - (self.rj*self.f2_func(X) - self.rm*self.f6_func(X))
                ) * np.exp(self.b_spr*(self.rm*X[4] - self.rj*X[0]))
            )
        )

    # def Lfhs(self,X):
    #     return(
    #         (self.rj**3)*self.k_spr*(self.b_spr**2)*(
    #             np.exp(self.b_spr*(self.rm*X[2] + self.rj*X[0]))
    #             - np.exp(self.b_spr*(self.rm*X[4] - self.rj*X[0]))
    #         ) * self.f1
    #         + (self.rj**2)*self.rm*self.k_spr*(self.b_spr**2)*(
    #             np.exp(self.b_spr*(self.rm*X[2] + self.rj*X[0]))
    #         ) * self.f3
    #         + (self.rj**2)*self.rm*self.k_spr*(self.b_spr**2)*(
    #             np.exp(self.b_spr*(self.rm*X[4] - self.rj*X[0]))
    #         ) * self.f5
    #     )
    # def Lf2hs(self,X):
    #     return(
    #         (
    #             self.rj**4*self.k_spr*self.b_spr**3*(
    #                 np.exp(
    #                     self.b_spr*(
    #                         self.rm*X[2]
    #                         + self.rj*X[0]
    #                     )
    #                 )
    #                 + np.exp(
    #                     self.b_spr*(
    #                         self.rm*X[4]
    #                         - self.rj*X[0]
    #                     )
    #                 )
    #             ) * self.f1(X)
    #             + self.rj**3*self.rm*self.k_spr*self.b_spr**3*(
    #                 np.exp(
    #                     self.b_spr*(
    #                         self.rm*X[2]
    #                         + self.rj*X[0]
    #                     )
    #                 )
    #             ) * self.f3(X)
    #             - self.rj**3*self.rm*self.k_spr*self.b_spr**3*(
    #                 np.exp(
    #                     self.b_spr*(
    #                         self.rm*X[4]
    #                         - self.rj*X[0]
    #                     )
    #                 )
    #             ) * self.f5(X)
    #         ) * self.f1(X)
    #         + (
    #             self.rj**3*self.k_spr*self.b_spr**2*(
    #                 np.exp(
    #                     self.b_spr*(
    #                         self.rm*X[2]
    #                         + self.rj*X[0]
    #                     )
    #                 )
    #                 - np.exp(
    #                     self.b_spr*(
    #                         self.rm*X[4]
    #                         - self.rj*X[0]
    #                     )
    #                 )
    #             )
    #         ) * self.f2(X)
    #         + (
    #             self.rj**3*self.rm*self.k_spr*self.b_spr**3*(
    #                 np.exp(
    #                     self.b_spr*(
    #                         self.rm*X[2]
    #                         + self.rj*X[0]
    #                     )
    #                 )
    #             ) * self.f1(X)
    #             + self.rj**2*self.rm**2*self.k_spr*self.b_spr**3*(
    #                 np.exp(
    #                     self.b_spr*(
    #                         self.rm*X[2]
    #                         + self.rj*X[0]
    #                     )
    #                 )
    #             ) * self.f3(X)
    #         ) * self.f3(X)
    #         + (
    #             self.rj**2*self.rm*self.k_spr*self.b_spr**2*(
    #                 np.exp(
    #                     self.b_spr*(
    #                         self.rm*X[2]
    #                         + self.rj*X[0]
    #                     )
    #                 )
    #             )
    #         ) * self.f4(X)
    #         + (
    #             self.rj**3*self.rm*self.k_spr*self.b_spr**3*(
    #                 -np.exp(
    #                     self.b_spr*(
    #                         self.rm*X[4]
    #                         - self.rj*X[0]
    #                     )
    #                 )
    #             ) * self.f1(X)
    #             + self.rj**2*self.rm**2*self.k_spr*self.b_spr**3*(
    #                 np.exp(
    #                     self.b_spr*(
    #                         self.rm*X[4]
    #                         - self.rj*X[0]
    #                     )
    #                 )
    #             ) * self.f5(X)
    #         ) * self.f5(X)
    #         + (
    #             self.rj**2*self.rm*self.k_spr*self.b_spr**2*(
    #                 np.exp(
    #                     self.b_spr*(
    #                         self.rm*X[4]
    #                         - self.rj*X[0]
    #                     )
    #                 )
    #             )
    #         ) * self.f6(X)
    #     )

    # def Phi(self,X):
    #     return(
    #         np.matrix([[
    #             self.h0(X),
    #             self.Lfh0(X),
    #             self.Lf2h0(X),
    #             self.Lf3h0(X),
    #             self.hs(X),
    #             self.Lfhs(X)
    #         ]]).T
    #     )
    def v0(self,X,x1d):
        k30 = 19.6
        k20 = 192.0
        k10 = 1101.9
        k00 = 3162.3
        assert np.all([k>0 for k in [k00,k10,k20,k30]]), "Error in gains for v0."
        assert k30*k20-k10 > 0, "Error in gains for v0."
        assert k30*k20*k10 - k30**2*k00 - k10 > 0, "Error in gains for v0."
        # import ipdb; ipdb.set_trace()
        result = (
            x1d[4]
            + k30*(x1d[3]-self.Lf3h0(X))
            + k20*(x1d[2]-self.Lf2h0(X))
            + k10*(x1d[1]-self.Lfh0(X))
            + k00*(x1d[0]-self.h0(X))
        )
        self.Lf4h0_list.append(result)
        return(result)
    def vs(self,X,Sd):
        k1s = 25.1
        k0s = 316.2
        result =(
            Sd[2]
            + k1s*(Sd[1]-self.Lfhs(X))
            + k0s*(Sd[0]-self.hs(X))
        )
        self.Lf2hs_list.append(result)
        return(result)

    def Q(self,X):
        B = np.matrix([
            [1/(self.Jm*self.Ij),0],
            [0,1/self.Jm]
        ])
        W = self.rj*self.rm*self.k_spr*self.b_spr*np.matrix([
            [
                np.exp(self.b_spr*(self.rm*X[2] + self.rj*X[0])),
                -np.exp(self.b_spr*(self.rm*X[4] - self.rj*X[0]))
            ],
            [
                self.rj*self.b_spr*(
                    np.exp(self.b_spr*(self.rm*X[2] + self.rj*X[0]))
                ),
                self.rj*self.b_spr*(
                    np.exp(self.b_spr*(self.rm*X[4] - self.rj*X[0]))
                )
            ]
        ])
        return(B*W)
    def return_input(self,X,x1d,Sd):
        """
            grad_f2_func = deriv_funcs["Gradient of f2"]
            hess_f2_func = deriv_funcs["Hessian of f2"]
            jac_f_func = deriv_funcs["Jacobian of f"]
        """
        try:
            # Q_inv = np.linalg.inv(self.Q(X))
            Q_inv = self.Q(X)**(-1)
        except:
            import ipdb; ipdb.set_trace()
        # import ipdb; ipdb.set_trace()
        return(
            Q_inv
            * (
                np.matrix([[-self.Lf4h0(X),-self.Lf2hs(X)]]).T
                + np.matrix([[self.v0(X,x1d),self.vs(X,Sd)]]).T
            )
        )

    def forward_simulation_FL(self,Time,X_o,X1d,Sd):
        assert len(X_o)==6, "X_o must have 6 elements, not " + str(len(X_o)) + "."
        dt = Time[1]-Time[0]
        U = np.zeros((2,len(Time)-1))
        X = np.zeros((6,len(Time)))
        Y = np.zeros((2,len(Time)))
        X[:,0] = X_o
        Y[:,0] = self.h(X[:,0])
        # self.update_state_variables(X_o)
        # deriv_funcs = self.create_derivative_functions()
        # self.grad_f2_func = deriv_funcs["Gradient of f2"]
        # self.hess_f2_func = deriv_funcs["Hessian of f2"]
        # self.jac_f_func = deriv_funcs["Jacobian of f"]
        statusbar=dsb(0,len(Time)-1,title="Forward Simulation (FBL)")
        for i in range(len(Time)-1):
            U[:,i] = (self.return_input(X[:,i],X1d[:,i],Sd[:,i])).flatten()
            # import ipdb; ipdb.set_trace()
            X[:,i+1] = (
                X[:,i]
                + dt*(
                    self.f(X[:,i])
                    + self.g(X[:,i])@U[:,np.newaxis,i]
                ).T
            )
            Y[:,i+1] = self.h(X[:,i+1])
            # self.update_state_variables(X[:,i+1])
            statusbar.update(i)
        return(X,U,Y)

def test_plant():
    plant1 = plant_pendulum_1DOF2DOF(**params)
    plant2 = plant_pendulum_1DOF2DOF(**params)
    dt = 0.001
    Time = np.arange(0,10+dt,dt)
    So = 100
    x1o = -np.pi/2
    x35o = (
        (1/(plant1.rm*plant1.b_spr))*np.log(
            So
            /(
                plant1.k_spr*plant1.b_spr*plant1.rj**2*(
                    np.exp(plant1.b_spr*plant1.rj*x1o)
                    + np.exp(-plant1.b_spr*plant1.rj*x1o)
                )
            )
        )
    )
    X_o = [x1o,0,x35o,0,x35o,0]
    X,U,Y = plant1.forward_simulation(Time,X_o)
    X1d = np.zeros((5,len(Time)))
    # X1d[0,int(2/dt):int(5/dt)] = -np.pi/2+np.ones((1,int((5-2)/dt)))
    # X1d[0,int(5/dt):int(8/dt)] = -np.pi/2-np.ones((1,int((5-2)/dt)))
    X1d[0,:] = -np.pi/4*np.ones((1,len(Time)))
    X1d[0,:] = LP_filt(1000, X1d[0,:])
    X1d[1,:] = np.gradient(X1d[0,:],dt)
    X1d[2,:] = np.gradient(X1d[1,:],dt)
    X1d[3,:] = np.gradient(X1d[2,:],dt)
    X1d[4,:] = np.gradient(X1d[3,:],dt)
    Sd = np.zeros((3,len(Time)))
    Sd[0,:] = So*np.ones((1,len(Time)))
    Sd[0,int(2/dt):int(5/dt)] = So+10*np.ones((1,int((5-2)/dt)))
    Sd[0,int(5/dt):int(8/dt)] = So-10*np.ones((1,int((5-2)/dt)))
    # Sd[0,:] = 20*np.cos(np.pi*Time)+(So-20)0.9*
    # Sd[1,:] = -20*np.pi*np.sin(np.pi*Time)0.9*
    X_FBL,U_FBL,Y_FBL = plant2.forward_simulation_FL(Time,X_o,X1d,Sd)
    plt.figure(figsize=(10,8))
    ax1=plt.gca()
    ax1.plot(Time,(180/np.pi)*Y[0,:].T,c="C0",linestyle=":")
    ax1.plot(Time,(180/np.pi)*Y_FBL[0,:].T,c="C0")
    ax1.plot(Time,(180/np.pi)*X1d[0,:],c="C0",linestyle="--")
    ax1.set_title(r"$-$ Actual; --- Desired; $\cdots$ Unforced", fontsize=16)
    ax1.set_xlabel("Time (s)")
    ax1.tick_params(axis='y', labelcolor="C0")
    ax1.set_ylabel('Position (deg.)', color="C0")
    ax2 = ax1.twinx()
    ax2.plot(Time,Y[1,:].T,c="C1",linestyle=":")
    ax2.plot(Time,Y_FBL[1,:].T,c="C1")
    ax2.plot(Time,Sd[0,:],c="C1",linestyle="--")
    ax2.tick_params(axis='y', labelcolor="C1")
    ax2.set_ylabel('Stiffness (Nm/rad.)', color="C1")

    plt.figure()
    plt.plot(Time,Y_FBL[1,:] - Sd[0,:])

    plt.show()
    return(Time,X_FBL,U_FBL,Y_FBL,plant1,plant2)
