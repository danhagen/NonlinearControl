from pendulum_eqns.integrator_backstepping_equations import *
from pendulum_eqns.initial_tension import *

if g == 0:
	MaxStep_MuscleVelocity = 0.05 # percentage of positive maximum.
	MuscleVelocity_Bounds =[[-5*lo1,5*lo1],[-5*lo2,5*lo2]]
else:
	MaxStep_MuscleVelocity = 1 # percentage of positive maximum.
	MuscleVelocity_Bounds =[[-5*lo1,5*lo1],[-1*lo2,1*lo2]]

def return_constraint_variables_muscle_velocity_driven(t,X):
	Coefficient1 = c2*c3*R1(X)*KT_1(X)
	Coefficient2 = c2*c4*R2(X)*KT_2(X)
	Constraint = A3(t,X,InitGlobals=True)
	return(Coefficient1,Coefficient2,Constraint)

def find_initial_values_Vm(**kwargs):
	X_o = np.array([Amp+Base,0])
	T = return_initial_tension(X_o,**kwargs)
	return(np.array([Amp+Base,0,T[0],T[1]]))

def animate_muscle_velocity_driven(t,X,U,**kwargs):
	"""
	Takes in Time (t - numpy.ndarray of shape (N,)), the state array (X - numpy.ndarray of shape (4,N)), and the input array (U - numpy.ndarray of shape (2,N)) and animates constraint equation over time.

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

	assert np.shape(X) == (4,len(t)) and str(type(X)) == "<class 'numpy.ndarray'>", "X must be a (4,N) numpy.ndarray"

	Bounds = kwargs.get("Bounds",MuscleVelocity_Bounds)
	assert type(Bounds) == list and np.shape(Bounds) == (2,2), "Bounds for Muscle Velocity Control must be a (2,2) list."
	assert Bounds[0][0]<Bounds[0][1],"Each set of bounds must be in ascending order."
	assert Bounds[1][0]<Bounds[1][1],"Each set of bounds must be in ascending order."

	dt = t[1]-t[0]
	fig = plt.figure(figsize=(10,8))
	ax1 = plt.gca()

	DescriptiveTitle = "Plotting Constraints vs. Time\nMuscle Velocity Driven"

	ax1.set_title(DescriptiveTitle,Fontsize=20,y=1)

	#Hyperbolic Constraint/Bounding Constraints
	Input1 = list(np.linspace(Bounds[0][0],Bounds[0][1],1000001))
	Input1.remove(0)
	Input1 = np.array(Input1)
	ax1.plot(Input1,lo1*lo2*0.001**2/Input1,'r',lw=2)
	Input2 = list(np.linspace(Bounds[1][0],Bounds[1][1],1000001))
	Input2.remove(0)
	Input2 = np.array(Input2)
	ax1.plot(lo1*lo2*0.001**2/Input2,Input2,'r',lw=2)

	ax1.plot([Bounds[0][0],Bounds[0][1]],[Bounds[1][0],Bounds[1][0]],'k--')
	ax1.plot([Bounds[0][0],Bounds[0][1]],[Bounds[1][1],Bounds[1][1]],'k--')
	ax1.plot([Bounds[0][0],Bounds[0][0]],[Bounds[1][0],Bounds[1][1]],'k--')
	ax1.plot([Bounds[0][1],Bounds[0][1]],[Bounds[1][0],Bounds[1][1]],'k--')

	Coefficient1,Coefficient2,Constraint1 = return_constraint_variables_muscle_velocity_driven(t[0],X[:,0])
	if abs(Coefficient1) <= 1e-7:
		LowerBound = Bounds[0][0]
		UpperBound = Bounds[0][1]
		if Constraint1/Coefficient2 > 0:
			LowerBound = Bounds[0][0]
			UpperBound = (lo1*(0.001)*lo2*(0.001))/(Constraint1/Coefficient2)
		else:
			LowerBound = (lo1*(0.001)*lo2*(0.001))/(Constraint1/Coefficient2)
			UpperBound = Bounds[0][1]
		FeasibleInput1 = (UpperBound-LowerBound)*np.random.rand(1000) + LowerBound
		FeasibleInput2 = np.array([Constraint1/Coefficient2]*1000)
	elif abs(Coefficient2) <= 1e-7:
		LowerBound = Constraint1/Coefficient1
		UpperBound = Constraint1/Coefficient1
		FeasibleInput1 = np.array([Constraint1/Coefficient1]*1000)
		if Constraint1/Coefficient1 < 0:
			LowerBound = (lo1*(0.001)*lo2*(0.001))/(Constraint1/Coefficient1)
			UpperBound = Bounds[1][1]
		else:
			LowerBound = Bounds[1][0]
			UpperBound = (lo1*(0.001)*lo2*(0.001))/(Constraint1/Coefficient1)
		FeasibleInput2 = (UpperBound-LowerBound)*np.random.rand(1000) + LowerBound
	elif np.sign(-Coefficient1) == np.sign(Coefficient2): # DIFFERENT SIGNS
		if (Constraint1 - Coefficient1*Bounds[0][0])/Coefficient2 < Bounds[1][0]:
			LowerBound = (Constraint1-Coefficient2*Bounds[1][0])/Coefficient1
		else:
			LowerBound = Bounds[0][0]

		if (Constraint1 - Coefficient1*Bounds[0][1])/Coefficient2 > Bounds[1][1]:
			UpperBound = (Constraint1-Coefficient2*Bounds[1][1])/Coefficient1
		else:
			UpperBound = Bounds[0][1]
		assert UpperBound>=LowerBound, "Error creating bounds - UpperBound should always be greater than LowerBound."
		HyperbolicBounds = np.sort([(Constraint1 - \
										np.sqrt(Constraint1**2 - 4*Coefficient1*Coefficient2*(lo1*0.001)*(lo2*0.001)))\
											/(2*Coefficient1), \
								 	(Constraint1 + \
										np.sqrt(Constraint1**2 - 4*Coefficient1*Coefficient2*(lo1*0.001)*(lo2*0.001)))\
											/(2*Coefficient1)])
		LowerBound = max([LowerBound,HyperbolicBounds[0]])
		UpperBound = min([UpperBound,HyperbolicBounds[1]])
		FeasibleInput1 = (UpperBound-LowerBound)*np.random.rand(1000) + LowerBound
		FeasibleInput2 = np.array([Constraint1/Coefficient2 - (Coefficient1/Coefficient2)*el \
								for el in FeasibleInput1])
	else: # np.sign(-Coefficient1) != np.sign(Coefficient2) SAME SIGN
		if (Constraint1 - Coefficient1*Bounds[0][0])/Coefficient2 > Bounds[1][1]:
			LowerBound = (Constraint1-Coefficient2*Bounds[1][1])/Coefficient1
		else:
			LowerBound = Bounds[0][0]

		if (Constraint1 - Coefficient1*Bounds[0][1])/Coefficient2 < Bounds[1][0]:
			UpperBound = (Constraint1-Coefficient2*Bounds[1][0])/Coefficient1
		else:
			UpperBound = Bounds[0][1]
		assert UpperBound>=LowerBound, "Error creating bounds - UpperBound should always be greater than LowerBound."
		if Constraint1**2 - 4*Coefficient1*Coefficient2*(lo1*0.001)*(lo2*0.001) > 0:
			HyperbolicBounds = np.sort([(Constraint1 - \
											np.sqrt(Constraint1**2 - 4*Coefficient1*Coefficient2*(lo1*0.001)*(lo2*0.001)))\
												/(2*Coefficient1), \
									 	(Constraint1 + \
											np.sqrt(Constraint1**2 - 4*Coefficient1*Coefficient2*(lo1*0.001)*(lo2*0.001)))\
												/(2*Coefficient1)])

			assert LowerBound < HyperbolicBounds[0] or HyperbolicBounds[1] < UpperBound, "No feasible solutions."

			FeasibleInput1 = []
			while len(FeasibleInput1)<1000:
				Random1 = (UpperBound-LowerBound)*np.random.rand() + LowerBound
				if Random1<HyperbolicBounds[0] or Random1>HyperbolicBounds[1]: FeasibleInput1.append(Random1)
			FeasinbleInput1 = np.array(FeasibleInput1)
			FeasibleInput2 = np.array([Constraint1/Coefficient2 - (Coefficient1/Coefficient2)*el \
									for el in FeasibleInput1])
		else:
			FeasibleInput1 = (UpperBound-LowerBound)*np.random.rand(1000) + LowerBound
			FeasibleInput2 = np.array([Constraint1/Coefficient2 - (Coefficient1/Coefficient2)*el \
									for el in FeasibleInput1])
	feasible = plt.Circle((U[:,0]),radius=MaxStep_MuscleVelocity,Color='b',alpha=0.5)
	ax1.add_patch(feasible)
	cline, = plt.plot(FeasibleInput1,FeasibleInput2,'b',lw=2)
	TimeText = plt.text(0.1,0.1,"t = " + str(t[0]),fontsize=16)
	chosenpoint, = plt.plot(U[:,0],c='k',marker='o')
	ax1.set_xlabel(r'$v_{m,1}$',fontsize=14)
	ax1.set_ylabel(r'$v_{m,2}$',fontsize=14)
	ax1.set_xlim([Bounds[0][0]-0.10*(np.diff(Bounds[0])[0]/2),\
					Bounds[0][1]+0.10*(np.diff(Bounds[0])[0]/2)])
	ax1.set_ylim([Bounds[1][0]-0.10*(np.diff(Bounds[1])[0]/2),\
					Bounds[1][1]+0.10*(np.diff(Bounds[1])[0]/2)])
	ax1.spines['right'].set_visible(False)
	ax1.spines['top'].set_visible(False)
	ax1.set_aspect('equal')

	def animate(i):
		Coefficient1,Coefficient2,Constraint1 = return_constraint_variables_muscle_velocity_driven(t[i],X[:,i])
		if abs(Coefficient1) <= 1e-7:
			LowerBound = Bounds[0][0]
			UpperBound = Bounds[0][1]
			if Constraint1/Coefficient2 > 0:
				LowerBound = Bounds[0][0]
				UpperBound = (lo1*(0.001)*lo2*(0.001))/(Constraint1/Coefficient2)
			else:
				LowerBound = (lo1*(0.001)*lo2*(0.001))/(Constraint1/Coefficient2)
				UpperBound = Bounds[0][1]
			FeasibleInput1 = (UpperBound-LowerBound)*np.random.rand(1000) + LowerBound
			FeasibleInput2 = np.array([Constraint1/Coefficient2]*1000)
		elif abs(Coefficient2) <= 1e-7:
			LowerBound = Constraint1/Coefficient1
			UpperBound = Constraint1/Coefficient1
			FeasibleInput1 = np.array([Constraint1/Coefficient1]*1000)
			if Constraint1/Coefficient1 < 0:
				LowerBound = (lo1*(0.001)*lo2*(0.001))/(Constraint1/Coefficient1)
				UpperBound = Bounds[1][1]
			else:
				LowerBound = Bounds[1][0]
				UpperBound = (lo1*(0.001)*lo2*(0.001))/(Constraint1/Coefficient1)
			FeasibleInput2 = (UpperBound-LowerBound)*np.random.rand(1000) + LowerBound
		elif np.sign(-Coefficient1) == np.sign(Coefficient2): # DIFFERENT SIGNS
			if (Constraint1 - Coefficient1*Bounds[0][0])/Coefficient2 < Bounds[1][0]:
				LowerBound = (Constraint1-Coefficient2*Bounds[1][0])/Coefficient1
			else:
				LowerBound = Bounds[0][0]

			if (Constraint1 - Coefficient1*Bounds[0][1])/Coefficient2 > Bounds[1][1]:
				UpperBound = (Constraint1-Coefficient2*Bounds[1][1])/Coefficient1
			else:
				UpperBound = Bounds[0][1]
			assert UpperBound>=LowerBound, "Error creating bounds - UpperBound should always be greater than LowerBound."
			HyperbolicBounds = np.sort([(Constraint1 - \
											np.sqrt(Constraint1**2 - 4*Coefficient1*Coefficient2*(lo1*0.001)*(lo2*0.001)))\
												/(2*Coefficient1), \
									 	(Constraint1 + \
											np.sqrt(Constraint1**2 - 4*Coefficient1*Coefficient2*(lo1*0.001)*(lo2*0.001)))\
												/(2*Coefficient1)])
			LowerBound = max([LowerBound,HyperbolicBounds[0]])
			UpperBound = min([UpperBound,HyperbolicBounds[1]])
			FeasibleInput1 = (UpperBound-LowerBound)*np.random.rand(1000) + LowerBound
			FeasibleInput2 = np.array([Constraint1/Coefficient2 - (Coefficient1/Coefficient2)*el \
									for el in FeasibleInput1])
		else: # np.sign(-Coefficient1) != np.sign(Coefficient2) SAME SIGN
			if (Constraint1 - Coefficient1*Bounds[0][0])/Coefficient2 > Bounds[1][1]:
				LowerBound = (Constraint1-Coefficient2*Bounds[1][1])/Coefficient1
			else:
				LowerBound = Bounds[0][0]

			if (Constraint1 - Coefficient1*Bounds[0][1])/Coefficient2 < Bounds[1][0]:
				UpperBound = (Constraint1-Coefficient2*Bounds[1][0])/Coefficient1
			else:
				UpperBound = Bounds[0][1]
			assert UpperBound>=LowerBound, "Error creating bounds - UpperBound should always be greater than LowerBound."
			if Constraint1**2 - 4*Coefficient1*Coefficient2*(lo1*0.001)*(lo2*0.001) > 0:
				HyperbolicBounds = np.sort([(Constraint1 - \
												np.sqrt(Constraint1**2 - 4*Coefficient1*Coefficient2*(lo1*0.001)*(lo2*0.001)))\
													/(2*Coefficient1), \
										 	(Constraint1 + \
												np.sqrt(Constraint1**2 - 4*Coefficient1*Coefficient2*(lo1*0.001)*(lo2*0.001)))\
													/(2*Coefficient1)])

				assert LowerBound < HyperbolicBounds[0] or HyperbolicBounds[1] < UpperBound, "No feasible solutions."

				FeasibleInput1 = []
				while len(FeasibleInput1)<1000:
					Random1 = (UpperBound-LowerBound)*np.random.rand() + LowerBound
					if Random1<HyperbolicBounds[0] or Random1>HyperbolicBounds[1]: FeasibleInput1.append(Random1)
				FeasinbleInput1 = np.array(FeasibleInput1)
				FeasibleInput2 = np.array([Constraint1/Coefficient2 - (Coefficient1/Coefficient2)*el \
										for el in FeasibleInput1])
			else:
				FeasibleInput1 = (UpperBound-LowerBound)*np.random.rand(1000) + LowerBound
				FeasibleInput2 = np.array([Constraint1/Coefficient2 - (Coefficient1/Coefficient2)*el \
										for el in FeasibleInput1])
		feasible.center = (U[:,i])
		if i<10:
			feasible.radius = 10*MaxStep_MuscleVelocity
		else:
			feasible.radius = MaxStep_MuscleVelocity
		cline.set_xdata(FeasibleInput1)
		cline.set_ydata(FeasibleInput2)
		chosenpoint.set_xdata(U[0,i])
		chosenpoint.set_ydata(U[1,i])
		TimeText.set_text("t = " + str(t[i]))
		return feasible,cline,chosenpoint,TimeText,


	# Init only required for blitting to give a clean slate.
	def init():
		ax1.plot(Input1,lo1*lo2*0.001**2/Input1,'r',lw=2)
		ax1.plot([Bounds[0][0],Bounds[0][1]],[Bounds[1][0],Bounds[1][0]],'k--')
		ax1.plot([Bounds[0][0],Bounds[0][1]],[Bounds[1][1],Bounds[1][1]],'k--')
		ax1.plot([Bounds[0][0],Bounds[0][0]],[Bounds[1][0],Bounds[1][1]],'k--')
		ax1.plot([Bounds[0][1],Bounds[0][1]],[Bounds[1][0],Bounds[1][1]],'k--')
		Coefficient1,Coefficient2,Constraint1 = return_constraint_variables_muscle_velocity_driven(t[0],X[:,0])
		if abs(Coefficient1) <= 1e-7:
			LowerBound = Bounds[0][0]
			UpperBound = Bounds[0][1]
			if Constraint1/Coefficient2 > 0:
				LowerBound = Bounds[0][0]
				UpperBound = (lo1*(0.001)*lo2*(0.001))/(Constraint1/Coefficient2)
			else:
				LowerBound = (lo1*(0.001)*lo2*(0.001))/(Constraint1/Coefficient2)
				UpperBound = Bounds[0][1]
			FeasibleInput1 = (UpperBound-LowerBound)*np.random.rand(1000) + LowerBound
			FeasibleInput2 = np.array([Constraint1/Coefficient2]*1000)
		elif abs(Coefficient2) <= 1e-7:
			LowerBound = Constraint1/Coefficient1
			UpperBound = Constraint1/Coefficient1
			FeasibleInput1 = np.array([Constraint1/Coefficient1]*1000)
			if Constraint1/Coefficient1 < 0:
				LowerBound = (lo1*(0.001)*lo2*(0.001))/(Constraint1/Coefficient1)
				UpperBound = Bounds[1][1]
			else:
				LowerBound = Bounds[1][0]
				UpperBound = (lo1*(0.001)*lo2*(0.001))/(Constraint1/Coefficient1)
			FeasibleInput2 = (UpperBound-LowerBound)*np.random.rand(1000) + LowerBound
		elif np.sign(-Coefficient1) == np.sign(Coefficient2): # DIFFERENT SIGNS
			if (Constraint1 - Coefficient1*Bounds[0][0])/Coefficient2 < Bounds[1][0]:
				LowerBound = (Constraint1-Coefficient2*Bounds[1][0])/Coefficient1
			else:
				LowerBound = Bounds[0][0]

			if (Constraint1 - Coefficient1*Bounds[0][1])/Coefficient2 > Bounds[1][1]:
				UpperBound = (Constraint1-Coefficient2*Bounds[1][1])/Coefficient1
			else:
				UpperBound = Bounds[0][1]
			assert UpperBound>LowerBound, "Error creating bounds - UpperBound should always be greater than LowerBound."
			HyperbolicBounds = np.sort([(Constraint1 - \
											np.sqrt(Constraint1**2 - 4*Coefficient1*Coefficient2*(lo1*0.001)*(lo2*0.001)))\
												/(2*Coefficient1), \
									 	(Constraint1 + \
											np.sqrt(Constraint1**2 - 4*Coefficient1*Coefficient2*(lo1*0.001)*(lo2*0.001)))\
												/(2*Coefficient1)])
			LowerBound = max([LowerBound,HyperbolicBounds[0]])
			UpperBound = min([UpperBound,HyperbolicBounds[1]])
			FeasibleInput1 = (UpperBound-LowerBound)*np.random.rand(1000) + LowerBound
			FeasibleInput2 = np.array([Constraint1/Coefficient2 - (Coefficient1/Coefficient2)*el \
									for el in FeasibleInput1])
		else: # np.sign(-Coefficient1) != np.sign(Coefficient2) SAME SIGN
			if (Constraint1 - Coefficient1*Bounds[0][0])/Coefficient2 > Bounds[1][1]:
				LowerBound = (Constraint1-Coefficient2*Bounds[1][1])/Coefficient1
			else:
				LowerBound = Bounds[0][0]

			if (Constraint1 - Coefficient1*Bounds[0][1])/Coefficient2 < Bounds[1][0]:
				UpperBound = (Constraint1-Coefficient2*Bounds[1][0])/Coefficient1
			else:
				UpperBound = Bounds[0][1]
			assert UpperBound>=LowerBound, "Error creating bounds - UpperBound should always be greater than LowerBound."
			if Constraint1**2 - 4*Coefficient1*Coefficient2*(lo1*0.001)*(lo2*0.001) > 0:
				HyperbolicBounds = np.sort([(Constraint1 - \
												np.sqrt(Constraint1**2 - 4*Coefficient1*Coefficient2*(lo1*0.001)*(lo2*0.001)))\
													/(2*Coefficient1), \
										 	(Constraint1 + \
												np.sqrt(Constraint1**2 - 4*Coefficient1*Coefficient2*(lo1*0.001)*(lo2*0.001)))\
													/(2*Coefficient1)])

				assert LowerBound < HyperbolicBounds[0] or HyperbolicBounds[1] < UpperBound, "No feasible solutions."

				FeasibleInput1 = []
				while len(FeasibleInput1)<1000:
					Random1 = (UpperBound-LowerBound)*np.random.rand() + LowerBound
					if Random1<HyperbolicBounds[0] or Random1>HyperbolicBounds[1]: FeasibleInput1.append(Random1)
				FeasinbleInput1 = np.array(FeasibleInput1)
				FeasibleInput2 = np.array([Constraint1/Coefficient2 - (Coefficient1/Coefficient2)*el \
										for el in FeasibleInput1])
			else:
				FeasibleInput1 = (UpperBound-LowerBound)*np.random.rand(1000) + LowerBound
				FeasibleInput2 = np.array([Constraint1/Coefficient2 - (Coefficient1/Coefficient2)*el \
										for el in FeasibleInput1])
		feasible = plt.Circle((U[:,0]),radius=MaxStep_MuscleVelocity,Color='b',alpha=0.5)
		feasible.set_visible(False)
		cline, = plt.plot(FeasibleInput1,FeasibleInput2,'b',lw=2)
		cline.set_visible(False)
		chosenpoint, = plt.plot(U[:,0],c='k',marker='o')
		chosenpoint.set_visible(False)
		TimeText = plt.text(0.75,0.75,"t = " + str(t[0]),fontsize=16)
		TimeText.set_visible(False)
		return feasible,cline,chosenpoint,TimeText,

	ani = animation.FuncAnimation(fig, animate, np.arange(1, np.shape(X)[1],1), init_func=init,interval=1, blit=False)
	plt.show()

def plot_individual_constraint_versus_time_muscle_velocity_driven(t,X,**kwargs):
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

	assert np.shape(X)[0] == 4 and str(type(X)) == "<class 'numpy.ndarray'>", "X must be a (4,N) numpy.ndarray"

	Return = kwargs.get("Return",False)
	assert type(Return) == bool, "Return must be either True or False."

	DescriptiveTitle = "Plotting Coefficients/Constraints vs. Time\nMuscle Velocity Driven\n$A\cdot u_{1} + B\cdot u_{2} = C$"
	fig, (ax1,ax2,ax3) = plt.subplots(1,3,figsize=(12,5))
	plt.subplots_adjust(wspace=0.4,top=0.8)
	plt.suptitle(DescriptiveTitle,Fontsize=20,y=0.95)

	A,B,C = [],[],[]
	for i in range(np.shape(X)[1]):
		Coefficient1,Coefficient2,Constraint1 = return_constraint_variables_muscle_velocity_driven(t[i],X[:,i])
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

def plot_individual_coefficient2_versus_time_muscle_velocity_driven(t,X,**kwargs):
	"""
	B = c2⋅c4⋅R2(X)⋅KT_2(X)

	Takes in Time (t - numpy.ndarray of shape (N,)), the state array (X - numpy.ndarray of shape (4,N)) and plots the 2nd Coefficient of the Constraint Equation over time as well as its components.

	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	**kwargs
	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	1) Return - must be a bool. Determines whether or not a function handle is returned. Default is False.

	"""
	import numpy as np
	import matplotlib.pyplot as plt

	assert np.shape(X)[0] == 4 and str(type(X)) == "<class 'numpy.ndarray'>", "X must be a (4,N) numpy.ndarray"

	Return = kwargs.get("Return",False)
	assert type(Return) == bool, "Return must be either True or False."

	fig, (ax1,ax2,ax3) = plt.subplots(1,3,figsize=(12,5))
	plt.subplots_adjust(top=0.8,hspace=0.4,bottom=0.1,left=0.1,right=0.975,wspace=0.4)
	DescriptiveTitle = "Plotting $2^{nd}$ Coefficient vs. Time\n$B=c_{2}c_{4}R_{2}(\\vec{x}(t))K_{T,2}(\\vec{x}(t))$"
	plt.suptitle(DescriptiveTitle,Fontsize=20,y=0.95)

	r2,kt_2,B = [],[],[]
	for i in range(np.shape(X)[1]):
		_,Coefficient2,_ = return_constraint_variables_muscle_velocity_driven(t[i],X[:,i])
		B.append(Coefficient2)
		r2.append(R2(X[:,i]))
		kt_2.append(KT_2(X[:,i]))

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

	ax3.plot(t[:np.shape(X)[1]],B,'b',lw=2)
	ax3.spines['right'].set_visible(False)
	ax3.spines['top'].set_visible(False)
	ax3.set_ylabel(r"$2^{nd}$ Coefficient")
	ax3.set_xticks(ax1.get_xticks())
	ax3.set_xticklabels([""]*len(ax1.get_xticks()))

	if Return == True:
		return(fig)
	else:
		plt.show()

def plot_individual_coefficient1_versus_time_muscle_velocity_driven(t,X,**kwargs):
    	"""
    	A = c2⋅c3⋅R1(X)⋅KT_1(X)

    	Takes in Time (t - numpy.ndarray of shape (N,)), the state array (X - numpy.ndarray of shape (4,N)) and plots the 1st Coefficient of the Constraint Equation over time as well as its components.

    	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    	**kwargs
    	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    	1) Return - must be a bool. Determines whether or not a function handle is returned. Default is False.

    	"""
    	import numpy as np
    	import matplotlib.pyplot as plt

    	assert np.shape(X)[0] == 4 and str(type(X)) == "<class 'numpy.ndarray'>", "X must be a (4,N) numpy.ndarray"

    	Return = kwargs.get("Return",False)
    	assert type(Return) == bool, "Return must be either True or False."

    	fig, (ax1,ax2,ax3) = plt.subplots(1,3,figsize=(12,5))
    	plt.subplots_adjust(top=0.8,hspace=0.4,bottom=0.1,left=0.075,right=0.975,wspace=0.4)
    	DescriptiveTitle = "Plotting $1^{st}$ Coefficient vs. Time\n$A=c_{2}c_{3}R_{1}(\\vec{x}(t))K_{T,1}(\\vec{x}(t))$"
    	plt.suptitle(DescriptiveTitle,Fontsize=20,y=0.95)

    	r1,kt_1,B = [],[],[]
    	for i in range(np.shape(X)[1]):
    		Coefficient1,_,_ = return_constraint_variables_muscle_velocity_driven(t[i],X[:,i])
    		B.append(Coefficient1)
    		r1.append(R1(X[:,i]))
    		kt_1.append(KT_1(X[:,i]))

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

    	ax3.plot(t[:np.shape(X)[1]],B,'r',lw=2)
    	ax3.spines['right'].set_visible(False)
    	ax3.spines['top'].set_visible(False)
    	ax3.set_ylabel(r"$1^{st}$ Coefficient")
    	ax3.set_xticks(ax1.get_xticks())
    	ax3.set_xticklabels([""]*len(ax1.get_xticks()))

    	if Return == True:
    		return(fig)
    	else:
    		plt.show()
