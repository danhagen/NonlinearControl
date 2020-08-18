import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from danpy.sb import dsb
from danpy.useful_functions import save_figures,is_number
from scipy import signal
import numdifftools as nd
import scipy as sp
from params import *
# from animate import *

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

        self.k0 = params.get(
            "Position Gains",
            {
                0 : 3162.3,
                1 : 1101.9,
                2 : 192.0,
                3 : 19.6
            }
        )
        self.ks = params.get(
            "Stiffness Gains",
            {
                0 : 316.2,
                1 : 25.1
            }
        )

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
                np.exp(self.b_spr*(self.rm*X[2]-self.rj*X[0]))
                - 1
            )
        )
    def tendon_2_FL_func(self,X):
        return(
            self.k_spr*(
                np.exp(self.b_spr*(self.rm*X[4]+self.rj*X[0]))
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
                - self.Lcm*self.mj*gr*np.sin(X[0]) # gravitational torque
                + self.rj*self.k_spr * (
                    np.exp(self.b_spr*(self.rm*X[2] - self.rj*X[0]))
                    - np.exp(self.b_spr*(self.rm*X[4] + self.rj*X[0]))
                ) # total coupling torque between motors and joint
            )/self.Ij
        )
    def df2dx1_func(self,X):
        result = (
            (
                -self.dCdx1(X) # Coriolis and centrifugal torques (zero)
                - self.Lcm*self.mj*gr*np.cos(X[0]) # gravitational torque
                - (self.rj**2)*self.k_spr*self.b_spr * (
                    np.exp(self.b_spr*(self.rm*X[2] - self.rj*X[0]))
                    + np.exp(self.b_spr*(self.rm*X[4] + self.rj*X[0]))
                ) # total coupling torque between motors and joint
            )/self.Ij
        )
        return(result)
    def d2f2dx12_func(self,X):
        return(
            (
                -self.d2Cdx12(X) # Coriolis and centrifugal torques (zero)
                + self.Lcm*self.mj*gr*np.sin(X[0]) # gravitational torque
                + (self.rj**3)*self.k_spr*(self.b_spr**2) * (
                    np.exp(self.b_spr*(self.rm*X[2] - self.rj*X[0]))
                    - np.exp(self.b_spr*(self.rm*X[4] + self.rj*X[0]))
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
                np.exp(self.b_spr*(self.rm*X[2] - self.rj*X[0]))
            ) / self.Ij
        )
    def d2f2dx1x5_func(self,X):
        """
        This is equivalently dSdb/Ij
        """
        return(
            -(self.rj**2)*self.rm*self.k_spr*(self.b_spr**2) * (
                np.exp(self.b_spr*(self.rm*X[4] + self.rj*X[0]))
            ) / self.Ij
        )
    def df2dx2_func(self,X):
        result = (
            (
                -self.dCdx2(X) # Coriolis and centrifugal torques (zero)
                - self.bj # damping torque
            )/self.Ij
        )
        return(result)
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
        result = (
            self.rj*self.rm*self.k_spr*self.b_spr * (
                np.exp(self.b_spr*(self.rm*X[2] - self.rj*X[0]))
            ) / self.Ij
        )
        return(result)
    def d2f2dx32_func(self,X):
        return(
            self.rj*(self.rm**2)*self.k_spr*(self.b_spr**2) * (
                np.exp(self.b_spr*(self.rm*X[2] - self.rj*X[0]))
            ) / self.Ij
        )
    def df2dx5_func(self,X):
        """
        Equivalently, this is the negative value of -Q_{12}/Ij
        """
        result = (
            -self.rj*self.rm*self.k_spr*self.b_spr * (
                np.exp(self.b_spr*(self.rm*X[4] + self.rj*X[0]))
            ) / self.Ij
        )
        return(result)
    def d2f2dx52_func(self,X):
        return(
            -self.rj*(self.rm**2)*self.k_spr*(self.b_spr**2) * (
                np.exp(self.b_spr*(self.rm*X[4] + self.rj*X[0]))
            ) / self.Ij
        )

    def f3_func(self,X):
        return(X[3])

    def f4_func(self,X):
        return(
            (
                -self.bm*X[3]
                - self.rm*self.k_spr*(
                    np.exp(self.b_spr*(self.rm*X[2] - self.rj*X[0]))
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
                    np.exp(self.b_spr*(self.rm*X[4] + self.rj*X[0]))
                    -1
                )
            )/self.Jm
        )

    def f(self,X):
        result = np.zeros((6,1))
        result[0,0] = self.f1
        result[1,0] = self.f2
        result[2,0] = self.f3
        result[3,0] = self.f4
        result[4,0] = self.f5
        result[5,0] = self.f6
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
            np.exp(self.b_spr*(self.rm*X[2] - self.rj*X[0]))
            + np.exp(self.b_spr*(self.rm*X[4] + self.rj*X[0]))
        )
        return(result)

    def forward_simulation(self,Time,X_o,U=None):
        """
        Building our own f_array to reduce the number of calls for f_funcs by making it a static call for each iteration in the FBL instance.
        """
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
        statusbar=dsb(0,len(Time)-1,title="Forward Simulation (Custom)")
        for i in range(len(Time)-1):
            f_array = np.zeros((6,1))
            f_array[0,0] = self.f1_func(X[:,i])
            f_array[1,0] = self.f2_func(X[:,i])
            f_array[2,0] = self.f3_func(X[:,i])
            f_array[3,0] = self.f4_func(X[:,i])
            f_array[4,0] = self.f5_func(X[:,i])
            f_array[5,0] = self.f6_func(X[:,i])
            X[:,i+1] = (
                X[:,i]
                + dt*(
                    f_array
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
        return(self.f2)
    def Lf3h0(self,X):
        result = (
            self.df2dx1*self.f1
            + self.df2dx2*self.f2
            + self.df2dx3*self.f3
            + self.df2dx5*self.f5
        )
        return(result)
    def Lf4h0(self,X):
        return(
            (
                self.d2f2dx12*self.f1
                + self.d2f2dx1x2*self.f2
                + self.df2dx2*self.df2dx1
                + self.d2f2dx1x3*self.f3
                + self.d2f2dx1x5*self.f5
            ) * self.f1
            + (
                self.d2f2dx1x2*self.f1
                + self.df2dx1
                + self.d2f2dx22*self.f2
                + (self.df2dx2**2)
            ) * self.f2
            + (
                self.d2f2dx1x3*self.f1
                + self.df2dx2*self.df2dx3
                + self.d2f2dx32*self.f3
            ) * self.f3
            + (
                self.df2dx3
            ) * self.f4
            + (
                self.d2f2dx1x5*self.f1
                + self.df2dx2*self.df2dx5
                + self.d2f2dx52*self.f5
            ) * self.f5
            + (
                self.df2dx5
            ) * self.f6
        )

    def hs(self,X):
        return(
            (self.rj**2)*self.k_spr*self.b_spr*(
                np.exp(self.b_spr*(self.rm*X[2] - self.rj*X[0]))
                + np.exp(self.b_spr*(self.rm*X[4] + self.rj*X[0]))
            )
        )
    def Lfhs(self,X):
        return(
            (self.rj**2)*self.k_spr*(self.b_spr**2)*(
                -(self.rj*self.f1 - self.rm*self.f3)*(
                    np.exp(self.b_spr*(self.rm*X[2] - self.rj*X[0]))
                )
                + (self.rj*self.f1 + self.rm*self.f5)*(
                    np.exp(self.b_spr*(self.rm*X[4] + self.rj*X[0]))
                )
            )
        )
    def Lf2hs(self,X):
        return(
            (self.rj**2)*self.k_spr*(self.b_spr**2)*(
                (
                    self.b_spr*(self.rj*self.f1 - self.rm*self.f3)**2
                    - self.rj*self.f2
                    + self.rm*self.f4
                ) * np.exp(self.b_spr*(self.rm*X[2] - self.rj*X[0]))
                + (
                    self.b_spr*(self.rj*self.f1 + self.rm*self.f5)**2
                    + self.rj*self.f2
                    + self.rm*self.f6
                ) * np.exp(self.b_spr*(self.rm*X[4] + self.rj*X[0]))
            )
        )

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
        result = (
            x1d[4]
            + self.k0[3]*(x1d[3]-self.Lf3h0(X))
            + self.k0[2]*(x1d[2]-self.Lf2h0(X))
            + self.k0[1]*(x1d[1]-self.Lfh0(X))
            + self.k0[0]*(x1d[0]-self.h0(X))
        )
        return(result)
    def vs(self,X,Sd):
        result =(
            Sd[2]
            + self.ks[1]*(Sd[1]-self.Lfhs(X))
            + self.ks[0]*(Sd[0]-self.hs(X))
        )
        return(result)

    def Q(self,X):
        B = np.matrix([
            [1/(self.Jm*self.Ij),0],
            [0,1/self.Jm]
        ])
        W = self.rj*self.rm*self.k_spr*self.b_spr*np.matrix([
            [
                np.exp(self.b_spr*(self.rm*X[2] - self.rj*X[0])),
                -np.exp(self.b_spr*(self.rm*X[4] + self.rj*X[0]))
            ],
            [
                self.rj*self.b_spr*(
                    np.exp(self.b_spr*(self.rm*X[2] - self.rj*X[0]))
                ),
                self.rj*self.b_spr*(
                    np.exp(self.b_spr*(self.rm*X[4] + self.rj*X[0]))
                )
            ]
        ])
        return(B*W)
    def return_input(self,X,x1d,Sd):
    	try:
    	    Q_inv = self.Q(X)**(-1)
    	except:
    	    import ipdb; ipdb.set_trace()
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
        U = np.zeros((2,len(Time)-1),dtype=np.float64)
        X = np.zeros((6,len(Time)),dtype=np.float64)
        X_measured = np.zeros((6,len(Time)),dtype=np.float64)
        Y = np.zeros((2,len(Time)),dtype=np.float64)
        X[:,0] = X_o
        Y[:,0] = self.h(X[:,0])
        self.update_state_variables(X_o)
        statusbar=dsb(0,len(Time)-1,title="Feedback Linearization")
        self.desiredOutput = np.array([X1d[0,:],Sd[0,:]])
        for i in range(len(Time)-1):
            if i>0:
                X_measured[0,i] = X[0,i]
                X_measured[1,i] = (X[0,i]-X[0,i-1])/self.dt
                X_measured[2,i] = X[2,i]
                X_measured[3,i] = (X[2,i]-X[2,i-1])/self.dt
                X_measured[4,i] = X[4,i]
                X_measured[5,i] = (X[4,i]-X[4,i-1])/self.dt
            else:
                X_measured[:,i] = X[:,i]
            U[:,i] = (self.return_input(X_measured[:,i],X1d[:,i],Sd[:,i])).flatten()
            X[:,i+1] = (
                X[:,i]
                + self.dt*(
                    self.f(X[:,i])
                    + self.g(X[:,i])@U[:,np.newaxis,i]
                ).T
            )
            Y[:,i+1] = self.h(X[:,i+1])
            self.update_state_variables(X[:,i+1])
            statusbar.update(i)
        return(X,U,Y,X_measured)

def test_plant():
    params["dt"]=0.001
    params["Simulation Duration"] = 100

    plant1 = plant_pendulum_1DOF2DOF(**params)
    plant2 = plant_pendulum_1DOF2DOF(**params)

    Time = np.arange(0,params["Simulation Duration"]+params["dt"],params["dt"])

    x1o = np.pi
    X_o = [x1o,0,plant2.rj*x1o/plant2.rm,0,-plant2.rj*x1o/plant2.rm,0]
    params["X_o"] = X_o

    X,U,Y = plant1.forward_simulation(Time,X_o)

    X1d = np.zeros((5,len(Time)))
    X1d[0,:] = np.pi*np.ones((1,len(Time)))
    # import ipdb; ipdb.set_trace()
    timeBreaks = [
        int(el*params["Simulation Duration"]/params["dt"])
        for el in [0, 0.13333, 0.21667, 0.41667, .57, .785, 1]
    ]

    X1d[0,timeBreaks[0]:timeBreaks[1]] =(
        np.pi*np.ones((1,int(np.diff(timeBreaks[0:2]))))
    )
    X1d[0,timeBreaks[1]:timeBreaks[2]] = (
        np.pi*np.ones((1,int(np.diff(timeBreaks[1:3]))))
        - 1
    )
    X1d[0,timeBreaks[2]:timeBreaks[3]] = (
        np.pi
        + 0.5*np.sin(3*np.pi*np.arange(0,20,params["dt"])/5)
    )
    X1d[0,timeBreaks[3]:timeBreaks[4]] = (
        np.pi*np.ones((1,int(np.diff(timeBreaks[3:5]))))
        + 1
    )
    X1d[0,timeBreaks[4]:timeBreaks[5]] = (
        np.pi*np.ones((1,int(np.diff(timeBreaks[4:6]))))
        + 0.5
    )
    X1d[0,timeBreaks[5]:timeBreaks[6]] =(
        np.pi*np.ones((1,int(np.diff(timeBreaks[5:]))))
    )
    X1d[0,:] = LP_filt(100, X1d[0,:])
    X1d[1,:] = np.gradient(X1d[0,:],params["dt"])
    X1d[2,:] = np.gradient(X1d[1,:],params["dt"])
    X1d[3,:] = np.gradient(X1d[2,:],params["dt"])
    X1d[4,:] = np.gradient(X1d[3,:],params["dt"])

    # Sd = np.zeros((3,len(Time)))
    # Sd[0,:] = 80 - 20*np.cos(16*np.pi*Time/25)
    # Sd[1,:] = 64*np.pi*np.sin(16*np.pi*Time/25)/5
    # Sd[2,:] = (4**5)*(np.pi**2)*np.cos(16*np.pi*Time/25)/(5**3)
    Sd = np.zeros((3,len(Time)))
    Sd[0,:] = 20*np.ones((1,len(Time)))


    X_FBL,U_FBL,Y_FBL,X_measured = plant2.forward_simulation_FL(Time,X_o,X1d,Sd)
    fig1 = plt.figure(figsize=(10,8))
    ax1=plt.gca()

    ax1.plot(Time,(180/np.pi)*Y_FBL[0,:].T,c="C0")
    ax1.plot(Time,(180/np.pi)*X1d[0,:],c="C0",linestyle="--")
    ax1.set_title(r"$-$ Actual; --- Desired", fontsize=16)
    ax1.set_xlabel("Time (s)")
    ax1.tick_params(axis='y', labelcolor="C0")
    ax1.set_ylabel('Position (deg.)', color="C0")
    # y1_min = np.floor((Y_FBL[0,:].min()*180/np.pi)/22.5)*22.5
    # y1_min = min([y1_min,np.floor((X1d[0,:].min()*180/np.pi)/22.5)*22.5])
    # y1_max = np.ceil((Y_FBL[0,:].max()*180/np.pi)/22.5)*22.5
    # y1_max = max([y1_max,np.ceil((X1d[0,:].max()*180/np.pi)/22.5)*22.5])
    y1_min = 0
    y1_max = 360
    yticks = np.arange(y1_min,y1_max+22.5,22.5)
    yticklabels = []
    for el in yticks:
    	if el%45==0:
    		yticklabels.append(str(int(el)) + r"$^\circ$")
    	else:
    		yticklabels.append("")
    ax1.set_yticks(yticks)
    ax1.set_yticklabels(yticklabels)
    ax2 = ax1.twinx()
    ax2.plot(Time,Y_FBL[1,:].T,c="C1")
    ax2.plot(Time,Sd[0,:],c="C1",linestyle="--")
    ax2.tick_params(axis='y', labelcolor="C1")
    ax2.set_ylabel('Stiffness (Nm/rad.)', color="C1")

    fig2 = plt.figure(figsize=(10,8))
    ax3=plt.gca()
    ax3.plot(Time,(180/np.pi)*(Y_FBL[0,:]-X1d[0,:]).T,c="C0")
    ax3.set_title("Error", fontsize=16)
    ax3.set_xlabel("Time (s)")
    ax3.tick_params(axis='y', labelcolor="C0")
    ax3.set_ylabel('Positional Error (deg.)', color="C0")
    yticklabels = [str(el)+r"$^\circ$" for el in ax3.get_yticks()]
    ax3.set_yticklabels(yticklabels)
    ax4 = ax3.twinx()
    ax4.plot(Time,Y_FBL[1,:] - Sd[0,:],c="C1")
    ax4.tick_params(axis='y', labelcolor="C1")
    ax4.set_ylabel('Stiffness Error (Nm/rad.)', color="C1")
    ax4.set_ylim([-0.1,0.1])
    ax4.set_yticks([-0.1,-0.05,0,0.05,0.1])

    fig3 = plt.figure(figsize=(10,8))
    tendonForce1Unforced = np.array(
        list(
            map(
                lambda X: plant1.tendon_1_FL_func(X),
                X.T
            )
        )
    )
    tendonDeformation1Unforced = np.array(
        list(
            map(
                lambda X: plant1.rm*X[2] - plant1.rj*X[0],
                X.T
            )
        )
    )
    tendonForce1_FBL = np.array(
        list(
            map(
                lambda X: plant2.tendon_1_FL_func(X),
                X_FBL.T
            )
        )
    )
    tendonDeformation1_FBL = np.array(
        list(
            map(
                lambda X: plant2.rm*X[2] - plant2.rj*X[0],
                X_FBL.T
            )
        )
    )

    tendonForce2Unforced = np.array(
        list(
            map(
                lambda X: plant1.tendon_2_FL_func(X),
                X.T
            )
        )
    )
    tendonDeformation2Unforced = np.array(
        list(
            map(
                lambda X: plant1.rm*X[4] + plant1.rj*X[0],
                X.T
            )
        )
    )
    tendonForce2_FBL = np.array(
        list(
            map(
                lambda X: plant2.tendon_2_FL_func(X),
                X_FBL.T
            )
        )
    )
    tendonDeformation2_FBL = np.array(
        list(
            map(
                lambda X: plant2.rm*X[4] + plant2.rj*X[0],
                X_FBL.T
            )
        )
    )

    minimumDeformation = min([
        tendonDeformation1_FBL.min(),
        tendonDeformation2_FBL.min()
    ])
    maximumDeformation = max([
        tendonDeformation1_FBL.max(),
        tendonDeformation2_FBL.max(),
        0.1
    ])
    deformationRange = maximumDeformation - minimumDeformation
    deformationArray = np.linspace(
        0,
        maximumDeformation+0.1*deformationRange,
        1001
    )
    actualForceLengthCurve = np.array(
        list(
            map(
                lambda x3: plant1.tendon_1_FL_func([0,0,x3/plant1.rm,0,0,0]),
                deformationArray
            )
        )
    )
    ax4 = fig3.add_subplot(211) # FL
    ax5 = fig3.add_subplot(212) # Time v Deformation
    ax4.plot(np.linspace(-1,0,1001),np.zeros((1001,)),'0.70')
    ax4.plot(deformationArray,actualForceLengthCurve,'0.70')
    ax4.plot(tendonDeformation1_FBL,tendonForce1_FBL,'r')
    ax4.plot(tendonDeformation2_FBL,tendonForce2_FBL,'g')
    ax4.set_xlim([
        minimumDeformation - 0.1*deformationRange,
        maximumDeformation + 0.1*deformationRange
    ])
    ax4.set_xlabel("Tendon Deformation (m)")
    ax4.set_ylabel("Tendon Tension (N)")
    ax4.spines['right'].set_visible(False)
    ax4.spines['top'].set_visible(False)

    ax5.plot(tendonDeformation1_FBL,-Time,'r')
    ax5.plot(tendonDeformation2_FBL,-Time,'g')
    ax5.set_ylabel("Time (s)")
    ax5.set_xlim([
        minimumDeformation - 0.1*deformationRange,
        maximumDeformation + 0.1*deformationRange
    ])
    ax4.set_xticklabels(["" for tick in ax4.get_xticks()])
    ax5.set_yticks([-Time[0],-Time[-1]])
    ax5.set_yticklabels([Time[0],Time[-1]])
    ax5.xaxis.tick_top()
    ax5.spines['right'].set_visible(False)
    ax5.spines['bottom'].set_visible(False)

    fig4 = plt.figure()
    ax6 = plt.gca()
    ax6.plot(Time,X_FBL[2,:]*180/np.pi,'r')
    ax6.plot(Time,X_FBL[4,:]*180/np.pi,'g')
    ax6.set_xlabel("Time (s)")
    ax6.set_ylabel("Motor Angles (deg)")
    ax6.spines["right"].set_visible(False)
    ax6.spines["top"].set_visible(False)

    # animate_trajectory(Time,X_FBL,U_FBL,Y_FBL,**params)
    plt.show()
    return(Time,X_FBL,U_FBL,Y_FBL,plant1,plant2)
