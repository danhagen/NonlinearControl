import numpy as np
import matplotlib.pyplot as plt
from danpy.sb import *
from matplotlib.patches import Ellipse
import matplotlib.animation as animation
import matplotlib.patches as patches
from params import *

def animate_trajectory(Time,X,U,**kwargs):
    SaveAsGif = kwargs.get("SaveAsGif",False)
    assert type(SaveAsGif)==bool, "SaveAsGif must be either True or False (Default)."

    FileName = kwargs.get("FileName","Pendulum_1DOF_2DOA_TT")
    assert type(FileName)==str,"FileName must be a str."

        # Angles must be in degrees for animation

    X1d = X[0,:]*(180/np.pi)
    X2d = X[1,:]*(180/np.pi)


    fig = plt.figure(figsize=(12,10))
    ax0 = plt.subplot2grid((2,4),(0,0),colspan=2) # animation
    ax1 = plt.subplot2grid((2,4),(0,2),colspan=2) # input
    ax2 = plt.subplot2grid((2,4),(1,0),colspan=2) # pendulum angle
    ax4 = plt.subplot2grid((2,4),(1,2),colspan=2) # pendulum angular velocity

    plt.suptitle("Cart-Pendulum Example",Fontsize=28,y=0.95)
    Pendulum_Width = 0.01*L1
    Pendulum_Length = L1

    Ground = plt.Rectangle(
                (-52*Pendulum_Width/4,-Pendulum_Length/4),
                52*Pendulum_Width/4,
                Pendulum_Length/2,
                Color='#4682b4')
    ax0.add_patch(Ground)


    Pendulum, = ax0.plot(
                    [
                        0,
                        Pendulum_Length*np.sin(X[0,0])
                    ],
                    [
                        0,
                        -Pendulum_Length*np.cos(X[0,0])
                    ],
                    Color='0.50',
                    lw = 10,
                    solid_capstyle='round'
                    )

    max_tau = U.max()
    if max_tau==0: max_tau=1

    k = 0.075*Pendulum_Length
    tau1_indicator, = ax0.plot(
                    0.75*Pendulum_Length*np.sin(
                            np.linspace(
                                1.05*X[0,0],
                                1.05*X[0,0] + (45*np.pi/180)*U[0,0]/max_tau,
                                20
                            )
                        ),
                    -0.75*Pendulum_Length*np.cos(
                            np.linspace(
                                1.05*X[0,0],
                                1.05*X[0,0] + (45*np.pi/180)*U[0,0]/max_tau,
                                20
                            )
                        ),
                    Color='r',
                    lw = 2,
                    solid_capstyle = 'round'
                    )
    tau1_indicator_arrow, = ax0.plot(
                    0.75*Pendulum_Length*np.sin(1.05*X[0,0] + (45*np.pi/180)*U[0,0]/max_tau)
                    + [
                        -k*np.sin((120*np.pi/180) - 1.05*X[0,0] - (45*np.pi/180)*U[0,0]/max_tau),
                        0,
                        -k*np.sin((60*np.pi/180) - 1.05*X[0,0] - (45*np.pi/180)*U[0,0]/max_tau)
                    ],
                    -0.75*Pendulum_Length*np.cos(1.05*X[0,0] + (45*np.pi/180)*U[0,0]/max_tau)
                    + [
                        -k*np.cos((120*np.pi/180) - 1.05*X[0,0] - (45*np.pi/180)*U[0,0]/max_tau),
                        0,
                        -k*np.cos((60*np.pi/180) - 1.05*X[0,0] - (45*np.pi/180)*U[0,0]/max_tau)
                    ],
                    Color='r',
                    lw = 2,
                    solid_capstyle='round'
                    )

    tau2_indicator, = ax0.plot(
                    0.75*Pendulum_Length*np.sin(
                            np.linspace(
                                0.95*X[0,0]-(45*np.pi/180)*U[1,0]/max_tau,
                                0.95*X[0,0],
                                20
                            )
                        ),
                    -0.75*Pendulum_Length*np.cos(
                            np.linspace(
                                0.95*X[0,0]-(45*np.pi/180)*U[1,0]/max_tau,
                                0.95*X[0,0],
                                20
                            )
                        ),
                    Color='g',
                    lw = 2,
                    solid_capstyle = 'round'
                    )
    tau2_indicator_arrow, = ax0.plot(
                    0.75*Pendulum_Length*np.sin(0.95*X[0,0] - (45*np.pi/180)*U[1,0]/max_tau)
                    + [
                        k*np.sin((60*np.pi/180) + 0.95*X[0,0] - (45*np.pi/180)*U[1,0]/max_tau),
                        0,
                        k*np.sin((120*np.pi/180) + 0.95*X[0,0] - (45*np.pi/180)*U[1,0]/max_tau)
                    ],
                    -0.75*Pendulum_Length*np.cos(0.95*X[0,0] - (45*np.pi/180)*U[1,0]/max_tau)
                    + [
                        -k*np.cos((60*np.pi/180) + 0.95*X[0,0] - (45*np.pi/180)*U[1,0]/max_tau),
                        0,
                        -k*np.cos((120*np.pi/180) + 0.95*X[0,0] - (45*np.pi/180)*U[1,0]/max_tau)
                    ],
                    Color='g',
                    lw = 2,
                    solid_capstyle='round'
                    )


    Pendulum_Attachment = plt.Circle((0,0),50*Pendulum_Width/4,Color='#4682b4')
    ax0.add_patch(Pendulum_Attachment)

    Pendulum_Rivet, = ax0.plot(
        [0],
        [0],
        c='k',
        marker='o',
        lw=2
        )

    ax0.get_xaxis().set_ticks([])
    ax0.get_yaxis().set_ticks([])
    ax0.set_frame_on(True)
    ax0.set_xlim([-1.10*Pendulum_Length,1.10*Pendulum_Length])
    ax0.set_ylim([-1.10*Pendulum_Length,1.10*Pendulum_Length])
    ax0.set_aspect('equal')


    TimeStamp = ax0.text(
        0,
        0.75*Pendulum_Length,
        "Time: "+str(Time[0])+" s",
        color='0.50',
        fontsize=16,
        horizontalalignment='center'
    )

    #Input

    Input1, = ax1.plot([0],[U[0,0]],color = 'r')
    Input2, = ax1.plot([0],[U[1,0]],color = 'g')
    ax1.set_xlim(0,Time[-1])
    ax1.set_xticks(list(np.linspace(0,Time[-1],5)))
    ax1.set_xticklabels([str(0),'','','',str(Time[-1])])
    if max(abs(U[0,:] - U[0,0]))<1e-7 and max(abs(U[1,:] - U[1,0]))<1e-7:
        ax1.set_ylim([min(U[:,0]) - 5,max(U[:,0]) + 5])
    else:
        RangeU = U.max()-U.min()
        ax1.set_ylim([U.min()-0.1*RangeU,U.max()+0.1*RangeU])

    ax1.spines['right'].set_visible(False)
    ax1.spines['top'].set_visible(False)
    ax1.set_title("Tendon Tensions (N)",fontsize=16,fontweight = 4,color = 'k',y = 0.95)

    #Pendulum Angle

    Angle, = ax2.plot([0],[X1d[0]],color = 'k')
    ax2.set_xlim(0,Time[-1])
    ax2.set_xticks(list(np.linspace(0,Time[-1],5)))
    ax2.set_xticklabels([str(0),'','','',str(Time[-1])])
    if max(abs(X1d-X1d[0]))<1e-7:
        ax2.set_ylim([X1d[0]-2,X1d[0]+2])
    else:
        RangeX1d= max(X1d)-min(X1d)
        ax2.set_ylim([min(X1d)-0.1*RangeX1d,max(X1d)+0.1*RangeX1d])
    ax2.spines['right'].set_visible(False)
    ax2.spines['top'].set_visible(False)
    ax2.set_title("Angle (deg)",fontsize=16,fontweight = 4,color = 'k',y = 0.95)

    # Angular Velocity

    AngularVelocity, = ax4.plot([0],[X2d[0]],color='k',linestyle='--')
    ax4.set_xlim(0,Time[-1])
    ax4.set_xticks(list(np.linspace(0,Time[-1],5)))
    ax4.set_xticklabels([str(0),'','','',str(Time[-1])])
    if max(abs(X2d-X2d[0]))<1e-7:
        ax4.set_ylim([X2d[0]-2,X2d[0]+2])
    else:
        RangeX2d= max(X2d)-min(X2d)
        ax4.set_ylim([min(X2d)-0.1*RangeX2d,max(X2d)+0.1*RangeX2d])
    ax4.spines['right'].set_visible(False)
    ax4.spines['top'].set_visible(False)
    ax4.set_title("Angular Velocity (deg/s)",fontsize=16,fontweight = 4,color = 'k',y = 0.95)

    def animate(i):
        Pendulum.set_xdata([0,Pendulum_Length*np.sin(X[0,i])])
        Pendulum.set_ydata([0,
                            -Pendulum_Length*np.cos(X[0,i])])
        tau1_indicator.set_xdata(
            0.75*Pendulum_Length*np.sin(
                np.linspace(
                    1.05*X[0,i],
                    1.05*X[0,i] + (45*np.pi/180)*U[0,i]/max_tau,
                    20
                )
            )
        )
        tau1_indicator.set_ydata(
            -0.75*Pendulum_Length*np.cos(
                np.linspace(
                    1.05*X[0,i],
                    1.05*X[0,i] + (45*np.pi/180)*U[0,i]/max_tau,
                    20
                )
            )
        )

        tau1_indicator_arrow.set_xdata(
            0.75*Pendulum_Length*np.sin(1.05*X[0,i] + (45*np.pi/180)*U[0,i]/max_tau)
            + [
                -k*np.sin((120*np.pi/180) - 1.05*X[0,i] - (45*np.pi/180)*U[0,i]/max_tau),
                0,
                -k*np.sin((60*np.pi/180) - 1.05*X[0,i] - (45*np.pi/180)*U[0,i]/max_tau)
            ]
        )
        tau1_indicator_arrow.set_ydata(
            -0.75*Pendulum_Length*np.cos(1.05*X[0,i] + (45*np.pi/180)*U[0,i]/max_tau)
            + [
                -k*np.cos((120*np.pi/180) - 1.05*X[0,i] - (45*np.pi/180)*U[0,i]/max_tau),
                0,
                -k*np.cos((60*np.pi/180) - 1.05*X[0,i] - (45*np.pi/180)*U[0,i]/max_tau)
            ]
        )

        tau2_indicator.set_xdata(
            0.75*Pendulum_Length*np.sin(
                np.linspace(
                    0.95*X[0,i]-(45*np.pi/180)*U[1,i]/max_tau,
                    0.95*X[0,i],
                    20
                )
            )
        )
        tau2_indicator.set_ydata(
            -0.75*Pendulum_Length*np.cos(
                np.linspace(
                    0.95*X[0,i]-(45*np.pi/180)*U[1,i]/max_tau,
                    0.95*X[0,i],
                    20
                )
            )
        )

        tau2_indicator_arrow.set_xdata(
            0.75*Pendulum_Length*np.sin(0.95*X[0,i] - (45*np.pi/180)*U[1,i]/max_tau)
            + [
                k*np.sin((60*np.pi/180) + 0.95*X[0,i] - (45*np.pi/180)*U[1,i]/max_tau),
                0,
                k*np.sin((120*np.pi/180) + 0.95*X[0,i] - (45*np.pi/180)*U[1,i]/max_tau)
            ]
        )
        tau2_indicator_arrow.set_ydata(
            -0.75*Pendulum_Length*np.cos(0.95*X[0,i] - (45*np.pi/180)*U[1,i]/max_tau)
            + [
                -k*np.cos((60*np.pi/180) + 0.95*X[0,i] - (45*np.pi/180)*U[1,i]/max_tau),
                0,
                -k*np.cos((120*np.pi/180) + 0.95*X[0,i] - (45*np.pi/180)*U[1,i]/max_tau)
            ]
        )


        TimeStamp.set_text("Time: "+"{:.2f}".format(Time[i])+" s",)

        Input1.set_xdata(Time[:i])
        Input1.set_ydata(U[0,:i])

        Input2.set_xdata(Time[:i])
        Input2.set_ydata(U[1,:i])

        Angle.set_xdata(Time[:i])
        Angle.set_ydata(X1d[:i])

        AngularVelocity.set_xdata(Time[:i])
        AngularVelocity.set_ydata(X2d[:i])

        return Pendulum,tau1_indicator,tau1_indicator_arrow,tau2_indicator,tau2_indicator_arrow,Input1,Input2,Angle,AngularVelocity,TimeStamp,

    # Init only required for blitting to give a clean slate.
    def init():
        Ground = plt.Rectangle(
                    (-52*Pendulum_Width/4,-Pendulum_Length/4),
                    52*Pendulum_Width/4,
                    Pendulum_Length/2,
                    Color='#4682b4')
        ax0.add_patch(Ground)


        Pendulum, = ax0.plot(
                        [
                            0,
                            Pendulum_Length*np.sin(X[0,0])
                        ],
                        [
                            0,
                            -Pendulum_Length*np.cos(X[0,0])
                        ],
                        Color='0.50',
                        lw = 10,
                        solid_capstyle='round'
                        )

        tau1_indicator, = ax0.plot(
                        0.75*Pendulum_Length*np.sin(
                                np.linspace(
                                    1.05*X[0,0],
                                    1.05*X[0,0] + (45*np.pi/180)*U[0,0]/max_tau,
                                    20
                                )
                            ),
                        -0.75*Pendulum_Length*np.cos(
                                np.linspace(
                                    1.05*X[0,0],
                                    1.05*X[0,0] + (45*np.pi/180)*U[0,0]/max_tau,
                                    20
                                )
                            ),
                        Color='r',
                        lw = 2,
                        solid_capstyle = 'round'
                        )
        tau1_indicator_arrow, = ax0.plot(
                        0.75*Pendulum_Length*np.sin(1.05*X[0,0] + (45*np.pi/180)*U[0,0]/max_tau)
                        + [
                            -k*np.sin((120*np.pi/180) - 1.05*X[0,0] - (45*np.pi/180)*U[0,0]/max_tau),
                            0,
                            -k*np.sin((60*np.pi/180) - 1.05*X[0,0] - (45*np.pi/180)*U[0,0]/max_tau)
                        ],
                        -0.75*Pendulum_Length*np.cos(1.05*X[0,0] + (45*np.pi/180)*U[0,0]/max_tau)
                        + [
                            -k*np.cos((120*np.pi/180) - 1.05*X[0,0] - (45*np.pi/180)*U[0,0]/max_tau),
                            0,
                            -k*np.cos((60*np.pi/180) - 1.05*X[0,0] - (45*np.pi/180)*U[0,0]/max_tau)
                        ],
                        Color='r',
                        lw = 2,
                        solid_capstyle='round'
                        )

        tau2_indicator, = ax0.plot(
                        0.75*Pendulum_Length*np.sin(
                                np.linspace(
                                    0.95*X[0,0]-(45*np.pi/180)*U[1,0]/max_tau,
                                    0.95*X[0,0],
                                    20
                                )
                            ),
                        -0.75*Pendulum_Length*np.cos(
                                np.linspace(
                                    0.95*X[0,0]-(45*np.pi/180)*U[1,0]/max_tau,
                                    0.95*X[0,0],
                                    20
                                )
                            ),
                        Color='g',
                        lw = 2,
                        solid_capstyle = 'round'
                        )
        tau2_indicator_arrow, = ax0.plot(
                        0.75*Pendulum_Length*np.sin(0.95*X[0,0] - (45*np.pi/180)*U[1,0]/max_tau)
                        + [
                            k*np.sin((60*np.pi/180) + 0.95*X[0,0] - (45*np.pi/180)*U[1,0]/max_tau),
                            0,
                            k*np.sin((120*np.pi/180) + 0.95*X[0,0] - (45*np.pi/180)*U[1,0]/max_tau)
                        ],
                        -0.75*Pendulum_Length*np.cos(0.95*X[0,0] - (45*np.pi/180)*U[1,0]/max_tau)
                        + [
                            -k*np.cos((60*np.pi/180) + 0.95*X[0,0] - (45*np.pi/180)*U[1,0]/max_tau),
                            0,
                            -k*np.cos((120*np.pi/180) + 0.95*X[0,0] - (45*np.pi/180)*U[1,0]/max_tau)
                        ],
                        Color='g',
                        lw = 2,
                        solid_capstyle='round'
                        )


        Pendulum_Attachment = plt.Circle((0,0),50*Pendulum_Width/4,Color='#4682b4')
        ax0.add_patch(Pendulum_Attachment)

        Pendulum_Rivet, = ax0.plot(
            [0],
            [0],
            c='k',
            marker='o',
            lw=2
            )

        TimeStamp = ax0.text(
            0,
            0.75*Pendulum_Length,
            "Time: "+"{:.2f}".format(Time[0])+" s",
            color='0.50',
            fontsize=16,
            horizontalalignment='center'
        )

        #Input

        Input1, = ax1.plot([0],[U[0,0]],color = 'r')
        Input2, = ax1.plot([0],[U[1,0]],color = 'g')

        #Pendulum Angle

        Angle, = ax2.plot([0],[X1d[0]],color = 'k')

        # Angular Velocity

        AngularVelocity, = ax4.plot([0],[X2d[0]],color='k',linestyle='--')


        Ground.set_visible(True)
        Pendulum.set_visible(False)
        tau1_indicator.set_visible(False)
        tau1_indicator_arrow.set_visible(False)
        tau2_indicator.set_visible(False)
        tau2_indicator_arrow.set_visible(False)
        Pendulum_Attachment.set_visible(False)
        Pendulum_Rivet.set_visible(False)
        TimeStamp.set_visible(False)
        Input1.set_visible(False)
        Input2.set_visible(False)
        Angle.set_visible(False)
        AngularVelocity.set_visible(False)

        return Ground,Pendulum,tau1_indicator,tau1_indicator_arrow,tau2_indicator,tau2_indicator_arrow,Pendulum_Attachment,Pendulum_Rivet,TimeStamp,Input1,Input2,Angle,AngularVelocity,

    dt = Time[1]-Time[0]
    if dt <= 0.0001:
        framestep=2000
    elif dt <= 0.001:
        framestep=200
    elif dt <= 0.01:
        framestep=10
    else:
        framestep=5
    ani = animation.FuncAnimation(fig, animate, frames=np.arange(0,len(Time)-1,framestep),init_func=init, blit=False)
    if SaveAsGif==True:
        ani.save(FileName+'.gif', writer='imagemagick', fps=10)
    plt.show()
