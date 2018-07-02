from integrator_backstepping_equations import *

if g == 0:
	MaxStep_Tension = 0.03 # percentage of positive maximum.
	Tension_Bounds = [[0,F_MAX1],[0,F_MAX2]]
else:
	MaxStep_Tension = 0.01 # percentage of positive maximum.
	Tension_Bounds = [[0,F_MAX1],[0,0.10*F_MAX2]]

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
