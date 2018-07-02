import numpy as np
import matplotlib.pyplot as plt
from termcolor import cprint,colored
from danpy.sb import dsb,get_terminal_width
from IB_muscle_velocities import *

def return_U_random_muscle_velocity(i,t:float,X,U,**kwargs):
	"""
	Takes in time scalar (float) (t), state numpy.ndarray (X) of shape (4,), and previous input numpy.ndarray (U) of shape (2,) and returns the input for this time step.

	Enforcing a hyperbolic domain constraint to allow for realistic lengthening/shortenting relationships.
	Input2 = (lo1*0.001)*(lo2*0.001)/Input1 = lo1*lo2/(10^6*Input1)

	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	**kwargs
	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	1) Noise - must be an numpy.ndarray of shape (2,). Default is np.zeros((1,2)).

	2) Seed - must be a scalar value. Default is None.

	3) Bounds - must be a (2,2) list with each row in ascending order. Default is given by MuscleVelocity_Bounds.

	4) MaxStep - must be a scalar (int or float). Default is MaxStep_MuscleVelocity.

	"""
	import random
	import numpy as np

	assert np.shape(X) == (4,) and str(type(X)) == "<class 'numpy.ndarray'>", "X must be a (4,) numpy.ndarray"
	assert np.shape(U) == (2,) and str(type(U)) == "<class 'numpy.ndarray'>", "U must be a (2,) numpy.ndarray"

	dt = t[1]-t[0]

	Noise = kwargs.get("Noise",np.zeros((2,)))
	assert np.shape(Noise) == (2,) and str(type(Noise)) == "<class 'numpy.ndarray'>", "Noise must be a (2,) numpy.ndarray"

	Seed = kwargs.get("Seed",None)
	assert type(Seed) in [float,int] or Seed is None, "Seed must be a float or an int or None."
	np.random.seed(Seed)

	Bounds = kwargs.get("Bounds",MuscleVelocity_Bounds)
	assert type(Bounds) == list and np.shape(Bounds) == (2,2), "Bounds for Muscle Velocity Control must be a (2,2) list."
	assert Bounds[0][0]<Bounds[0][1],"Each set of bounds must be in ascending order."
	assert Bounds[1][0]<Bounds[1][1],"Each set of bounds must be in ascending order."

	MaxStep = kwargs.get("MaxStep",MaxStep_MuscleVelocity)
	assert type(MaxStep) in [int,float], "MaxStep for Muscle Velocity Controller should be an int or float."

	Coefficient1,Coefficient2,Constraint1 = return_constraint_variables_muscle_velocity_driven(t[i],X)

	if Constraint1 != 0:
		assert Coefficient1!=0 or Coefficient2!=0, "Error with Coefficients. Shouldn't be zero with nonzero constraint."
	else:
		assert Coefficient1!=0 or Coefficient2!=0, "Error with Constraint. 0 = 0 implies all inputs valid."

	Roots = np.sort(\
				np.array(\
	     			list(\
	     				set(\
		     				np.roots(\
					     				[1,\
						     				-Constraint1/Coefficient1,\
						     					Coefficient2*lo1*lo2*(10**-6)/Coefficient1]\
																							)))))
	Roots = Roots[np.isreal(Roots)]

	if Coefficient1 == 0:
		if Constraint1/Coefficient2 > 0:
			LowerBound = Bounds[0][0]
			UpperBound = (lo1*(0.001)*lo2*(0.001))/(Constraint1/Coefficient2)
		else:
			LowerBound = (lo1*(0.001)*lo2*(0.001))/(Constraint1/Coefficient2)
			UpperBound = Bounds[0][1]
		FeasibleInput1 = (UpperBound-LowerBound)*np.random.rand(1000) + LowerBound
		FeasibleInput2 = np.array([Constraint1/Coefficient2]*1000)
	elif Coefficient2 == 0:
		FeasibleInput1 = np.array([Constraint1/Coefficient1]*1000)
		if Constraint1/Coefficient1 < 0:
			LowerBound = (lo1*(0.001)*lo2*(0.001))/(Constraint1/Coefficient1)
			UpperBound = Bounds[1][1]
		else:
			LowerBound = Bounds[1][0]
			UpperBound = (lo1*(0.001)*lo2*(0.001))/(Constraint1/Coefficient1)
		FeasibleInput2 = (UpperBound-LowerBound)*np.random.rand(1000) + LowerBound
	else:
		assert 0 not in Roots, "Zero should not be a root. (Implies Coefficient2 == 0)"
		if len(Roots) in [0,1]:
			SortedBounds = np.sort([(Constraint1-Coefficient2*Bounds[1][0])/Coefficient1,\
										(Constraint1-Coefficient2*Bounds[1][1])/Coefficient1])
			LowerBound = max(Bounds[0][0], SortedBounds[0])
			UpperBound = min(Bounds[0][1], SortedBounds[1])
			assert UpperBound >= LowerBound, "Error generating bounds. Not feasible!"
			FeasibleInput1 = (UpperBound-LowerBound)*np.random.rand(1000) + LowerBound
			FeasibleInput2 = np.array([Constraint1/Coefficient2 - (Coefficient1/Coefficient2)*el \
									for el in FeasibleInput1])
		elif (Roots<0).all() or (Roots>0).all():
			SortedBounds = np.sort([(Constraint1-Coefficient2*Bounds[1][0])/Coefficient1,\
										(Constraint1-Coefficient2*Bounds[1][1])/Coefficient1])
			LowerBound = max(Bounds[0][0], SortedBounds[0])
			UpperBound = min(Bounds[0][1], SortedBounds[1])
			ConstraintLength1 = Coefficient1/(2*Coefficient2)*(LowerBound**2-Roots[0]**2) \
									- Constraint1/Coefficient2*(LowerBound-Roots[0])
			ConstraintLength1 = ConstraintLength1*(ConstraintLength1>0)
			ConstraintLength2 = Coefficient1/(2*Coefficient2)*(Roots[1]**2-UpperBound**2) \
									- Constraint1/Coefficient2*(Roots[1]-UpperBound)
			ConstraintLength2 = ConstraintLength2*(ConstraintLength2>0)
			assert ConstraintLength1!=0 or ConstraintLength2!=0, \
								"Error generating bounds. Not feasible!"
			N1 = int(np.round(1000*ConstraintLength1/(ConstraintLength1+ConstraintLength2)))
			N2 = 1000-N1
			FeasibleInput1_1 = (Roots[0]-LowerBound)*np.random.rand(N1) + LowerBound
			FeasibleInput1_2 = (UpperBound-Roots[1])*np.random.rand(N2) + Roots[1]
			FeasibleInput1 = np.concatenate([FeasibleInput1_1,FeasibleInput1_2])
			FeasibleInput2 = np.array([Constraint1/Coefficient2 - (Coefficient1/Coefficient2)*el \
									for el in FeasibleInput1])
		else: # not((Roots<0).all()) and not((Roots>0).all()):
			SortedBounds = np.sort([(Constraint1-Coefficient2*Bounds[1][0])/Coefficient1,\
					(Constraint1-Coefficient2*Bounds[1][1])/Coefficient1])
			LowerBound = max(Bounds[0][0], SortedBounds[0],Roots[0])
			UpperBound = min(Bounds[0][1], SortedBounds[1],Roots[1])
			assert UpperBound >= LowerBound, "Error with Bounds. Infeasible!"
			FeasibleInput1 = (UpperBound-LowerBound)*np.random.rand(1000) + LowerBound
			FeasibleInput2 = np.array([Constraint1/Coefficient2 - (Coefficient1/Coefficient2)*el \
									for el in FeasibleInput1])

	def plot_constraints():
		import matplotlib.pyplot as plt
		plt.figure()
		Input1 = np.linspace(LowerBound,UpperBound,1001)
		Input2 = np.array([Constraint1/Coefficient2 - (Coefficient1/Coefficient2)*el \
								for el in Input1])
		plt.plot(Input1,Input2,'k--')
		plt.plot(Input1, (lo1*0.001)*(lo2*0.001)/Input1,'r')
		plt.scatter(FeasibleInput1,FeasibleInput2,c='g',marker = '.')
		plt.ylim(Bounds[1])
		plt.show()

	"""
	Checking to see which inputs have the appropriate allowable step size. In normalized muscle velocity.
	"""
	euclid_dist = np.array(list(map(lambda u1,u2: np.sqrt(((U[0]-u1)/lo1)**2 + ((U[1]-u2)/lo2)**2),\
							FeasibleInput1,FeasibleInput2)))
    import ipdb; ipdb.set_trace()
	feasible_index = np.where(euclid_dist <= \
									(MaxStep*(t[i]>=100*dt) + \
									 	10.0*MaxStep*(50*dt<=t[i]<100*dt) + \
											50.0*MaxStep*(t[i]<50*dt)))
	next_index = np.random.choice(feasible_index[0])
	u1 = FeasibleInput1[next_index]
	u2 = FeasibleInput2[next_index]
	return(np.array([u1,u2]))

def return_U_random_activations(i,t,X,U,**kwargs):
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

	MaxStep = kwargs.get("MaxStep",MaxStep_Activation)
	assert type(MaxStep) in [int,float], "MaxStep for Muscle Activation Controller should be an int or float."

	Coefficient1,Coefficient2,Constraint1 =\
	 			return_constraint_variables_muscle_activation_driven(t[i],X)
	assert Coefficient1!=0 and Coefficient2!=0, "Error with Coefficients. Shouldn't both be zero."
	if Constraint1 < 0:
		assert not(Coefficient1 > 0 and Coefficient2 > 0), "Infeasible activations. (Constraint1 < 0, Coefficient1 > 0, Coefficient2 > 0)"
	if Constraint1 > 0:
		assert not(Coefficient1 < 0 and Coefficient2 < 0), "Infeasible activations. (Constraint1 > 0, Coefficient1 < 0, Coefficient2 < 0)"

	AllowableBounds_x = np.array([U[0]-MaxStep,U[0]+MaxStep])
	AllowableBounds_y = np.array([U[1]-MaxStep,U[1]+MaxStep])

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
		# if UpperBound_x < LowerBound_x: import ipdb; ipdb.set_trace()
		assert UpperBound_x >= LowerBound_x, "Error generating bounds. Not feasible!"
		FeasibleInput1 = (UpperBound_x-LowerBound_x)*np.random.rand(1000) + LowerBound_x
		FeasibleInput2 = np.array([Constraint1/Coefficient2 - (Coefficient1/Coefficient2)*el \
								for el in FeasibleInput1])
	"""
	Checking to see which inputs have the appropriate allowable step size.
	"""

	next_index = np.random.choice(range(1000))
	u1 = FeasibleInput1[next_index]
	u2 = FeasibleInput2[next_index]
	return(np.array([u1,u2]))

def run_sim_rand_Vm(**kwargs):
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
        X = np.zeros((4,N))
        X[:,0] = find_initial_values_Vm()
        U = np.zeros((2,N))
        U[:,0] = [0,0]

        AddNoise = False
        if AddNoise == True:
            np.random.seed(seed=None)
            NoiseArray = np.random.normal(loc=0.0,scale=0.2,size=(2,N))
        else:
            NoiseArray = np.zeros((2,N))

        try:
            cprint("Attempt #" + str(int(AttemptNumber)) + ":\n", 'green')
            statusbar = dsb(0,N-1,title=run_sim_rand_Vm.__name__)
            for i in range(N-1):
                U[:,i+1] = return_U_random_muscle_velocity(i,Time,X[:,i],U[:,i],Noise = NoiseArray[:,i])
                X[:,i+1] = X[:,i] + dt*np.array([	dX1_dt(X[:,i]),\
                									dX2_dt(X[:,i]),\
                									dX3_dt(X[:,i],U=U[:,i+1]),\
                									dX4_dt(X[:,i],U=U[:,i+1])])
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

def run_N_sim_rand_Vm(**kwargs):
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
		TotalX[j],TotalU[j] = run_sim_rand_Vm(**kwargs)

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

def plot_N_sim_rand_Vm(t,TotalX,TotalU,**kwargs):
	Return = kwargs.get("Return",False)
	assert type(Return) == bool, "Return should either be True or False"

	fig1 = plt.figure(figsize = (9,7))
	fig1_title = "Underdetermined Forced-Pendulum Example"
	plt.title(fig1_title,fontsize=16,color='gray')
	statusbar = dsb(0,np.shape(TotalX)[0],title=(plot_N_sim_rand_Vm.__name__ + " (" + fig1_title +")"))
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
	statusbar.reset(title=(plot_N_sim_rand_Vm.__name__ + " (" + fig2_title +")"))
	for j in range(np.shape(TotalX)[0]):
		plt.plot(t, (r(t)-TotalX[j,0,:])*180/np.pi,color='0.70')
		statusbar.update(j)
	plt.xlabel("Time (s)")
	plt.ylabel("Error (Deg)")

	statusbar.reset(
		title=(
			plot_N_sim_rand_Vm.__name__
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
