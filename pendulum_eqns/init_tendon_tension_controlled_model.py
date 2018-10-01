from pendulum_eqns.integrator_backstepping_equations import *
from pendulum_eqns.initial_tension import *

def return_constraint_variables(t,X):
    Coefficient1 = c2*R1(X)
    Coefficient2 = c2*R2(X)
    Constraint = A2(t,X)
    return(Coefficient1,Coefficient2,Constraint)

def find_initial_values_TT(**kwargs):
	"""
	This function returns initial conditions for the system that starts from rest. (Ex. pendulum_eqns.reference_trajectories._01)
    """

	Seed = kwargs.get("Seed",None)
	assert type(Seed) in [float,int] or Seed is None, "Seed must be a float or an int or None."
	np.random.seed(Seed)

	X_o = np.array([Amp+Base,0])
	T = return_initial_tension(X_o,Bounds=[[0,0.05*F_MAX1],[0,0.05*F_MAX2]],InitialAngularAcceleration=d2r(0))
	return(X_o,T)

def animate_input_vs_time(t,X,U,**kwargs):
	"""
	Takes in Time (t - numpy.ndarray of shape (N,)), the state array (X - numpy.ndarray of shape (2,N)), and the input array (U - numpy.ndarray of shape (2,N)) and animates constraint equation over time.

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

	MaxStep = kwargs.get("MaxStep",MaxStep_Tension)
	assert type(MaxStep) in [int,float], "MaxStep for Muscle Activation Controller should be an int or float."

	Bounds = kwargs.get("Bounds",Tension_Bounds)
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

	Coefficient1,Coefficient2,Constraint1 = return_constraint_variables(t[0],X[:,0])

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
		Coefficient1,Coefficient2,Constraint1 = return_constraint_variables(t[i],X[:,i])
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
		Coefficient1,Coefficient2,Constraint1 = return_constraint_variables(t[0],X[:,0])
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
