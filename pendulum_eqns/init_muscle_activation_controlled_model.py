from pendulum_eqns.integrator_backstepping_equations import *
from pendulum_eqns.initial_tension import *
from pendulum_eqns.physiology.muscle_params_BIC_TRI import *

MaxStep_Activation = 1 # percentage of positive maximum (1)
Activation_Bounds = [[0,1],[0,1]]

def return_constraint_variables(t,X):
	Coefficient1 = c2*c3*c6*R1(X)*KT_1(X)*FLV_1(X)
	Coefficient2 = c2*c4*c10*R2(X)*KT_2(X)*FLV_2(X)
	Constraint = A4(t,X,InitGlobals=True)
	return(Coefficient1,Coefficient2,Constraint)

def return_random_initial_muscle_lengths_and_activations(InitialTension,**kwargs):
	PlotBool = kwargs.get("PlotBool",False)
	assert type(PlotBool)==bool,"PlotBool must be a boolean. Default is False."

	# L1 = np.linspace(0.5*lo1, 1.5*lo1, 1001)
	mu1, sigma1 = lo1, 0.1*lo1
	L1 = np.array(list(sorted(np.random.normal(mu1, sigma1, 1001))))
	U1 = np.array(
			list(
				map(
					lambda l:
						(InitialTension[0][0]/(F_MAX1*np.cos(α1))
						- F_PE1_1([0,0,0,0,l,0,0,0])
						)
						/ FL(l,lo1)
					,L1
				)
			)
		)
	thresh1 = (
			lo1
			* L_CE_max_1
			* (
				k_1
				* np.log(
					np.exp(
						InitialTension[0][0]
						/ (F_MAX1*c_1*k_1*np.cos(α1)))
						- 1
						)
					+ Lr1
				)
			)

	# L2 = np.linspace(0.5*lo2, 1.5*lo2, 1001)
	mu2, sigma2 = lo2, 0.1*lo2
	L2 = np.array(list(sorted(np.random.normal(mu2, sigma2, 1001))))
	U2 = np.array(
			list(
				map(
					lambda l:
						(InitialTension[1][0]/(F_MAX2*np.cos(α2))
						- F_PE1_2([0,0,0,0,0,l,0,0])
						)
						/ FL(l,lo2)
					,L2
				)
			)
		)
	thresh2 = (
			lo2
			* L_CE_max_2
			* (
				k_1
				* np.log(
					np.exp(
						InitialTension[1][0]
						/ (F_MAX2*c_1*k_1*np.cos(α2)))
						- 1
						)
					+ Lr1
				)
			)
	PositiveIndeces = np.intersect1d(
						np.where(L1<=thresh1),
						np.where(L2<=thresh2)
						)
	assert len(PositiveIndeces)>0, \
			("Error finding positive activations for given tension levels\nMuscle 1: "
			+ str(InitialTension[0][0])
			+ "\nMuscle 2: "
			+ str(InitialTension[1][0])
			)

	if PlotBool == True:
		plt.figure(figsize=(10,8))
		plt.title(r"Viable Initial $l_{m,1}$ and $u_{1}$ Values")
		plt.xlabel(r"$l_{m,1}$ (m)",fontsize=14)
		plt.ylabel(r"$u_{1}$",fontsize=14)
		L1_range = L1.max()-L1.min()
		L1_linspace = np.linspace(0.5*lo1,1.5*lo1,1001)
		U1_linspace = np.array(
				list(
					map(
						lambda l:
							(InitialTension[0][0]/(F_MAX1*np.cos(α1))
							- F_PE1_1([0,0,0,0,l,0,0,0])
							)
							/ FL(l,lo1)
						,L1_linspace
					)
				)
			)
		plt.plot(L1_linspace,U1_linspace,'k')
		plt.scatter(L1[PositiveIndeces],U1[PositiveIndeces])
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
		L2_range = L2.max()-L2.min()
		L2_linspace = np.linspace(0.5*lo2,1.5*lo2,1001)
		U2_linspace = np.array(
				list(
					map(
						lambda l:
							(InitialTension[1][0]/(F_MAX2*np.cos(α2))
							- F_PE1_2([0,0,0,0,0,l,0,0])
							)
							/ FL(l,lo2)
						,L2_linspace
					)
				)
			)
		plt.plot(L2_linspace,U2_linspace,'k')
		plt.scatter(L2[PositiveIndeces],U2[PositiveIndeces])
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
	return( L1[PositiveIndeces],
			U1[PositiveIndeces],
			L2[PositiveIndeces],
			U2[PositiveIndeces]
			)

def find_viable_initial_values(**kwargs):
	"""
	This function returns initial conditions for the system that starts from rest. (Ex. pendulum_eqns.reference_trajectories._01)

	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	**kwargs

	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	1) FixedInitialTension - Must be a (2,1) numpy.ndarray. Run find_initial_tension outside of the loop for a given seed and then feed it through the pipeline.

	2) ReturnAll - Can return all initial values for a given tension level. Will be fed through to return_random_initial_muscle_lengths_and_activations.

	3) Seed - Can see the random tension generated. When FixedInitialTension is provided, this seed will apply only to the initial conditions for activation and muscle length.
	"""
	FixedInitialTension = kwargs.get("FixedInitialTension",None)
	assert (FixedInitialTension is None) or \
			(str(type(FixedInitialTension)) == "<class 'numpy.ndarray'>"
			and np.shape(FixedInitialTension) == (2,1)),\
		(
		"FixedInitialTension must either be None (Default) or a (2,1) numpy.ndarray."
		+ "\nCurrent type: "
		+ str(type(FixedInitialTension))
		+ "\nCurrent shape: "
		+ str(np.shape(FixedInitialTension))
		)

	ReturnAll = kwargs.get("ReturnAll",False)
	assert type(ReturnAll)==bool, "ReturnAll must be a bool."

	Seed = kwargs.get("Seed",None)
	assert type(Seed) in [float,int] or Seed is None, "Seed must be a float or an int or None."
	np.random.seed(Seed)

	X_o = np.array([Amp+Base,0])
	if FixedInitialTension is None:
		T = return_initial_tension(X_o)
	else:
		T = FixedInitialTension
	L1,U1,L2,U2 = return_random_initial_muscle_lengths_and_activations(T,**kwargs)
	rand_index = np.random.choice(len(L1),2)

	if ReturnAll == False:
		return(
			T,
			np.array([L1[rand_index[0]],L2[rand_index[1]]]),
			np.array([U1[rand_index[0]],U2[rand_index[1]]])
			)
	else:
		return(T,L1,L2,U1,U2)

def animate_input_vs_time(t,X,U,**kwargs):
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

def plot_individual_constraint_versus_time(
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
		Coefficient1,Coefficient2,Constraint1 = return_constraint_variables(t[i],X[:,i])
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

def plot_individual_coefficient2_versus_time(
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
		_,Coefficient2,_ = return_constraint_variables(t[i],X[:,i])
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

def plot_individual_coefficient1_versus_time(
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
		Coefficient1,_,_ = return_constraint_variables(t[i],X[:,i])
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
