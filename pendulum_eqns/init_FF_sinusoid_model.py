from pendulum_eqns.integrator_backstepping_equations import *
from pendulum_eqns.initial_tension import *

ActivationBounds = [[0,1],[0,1]]

def return_random_initial_muscle_lengths_and_activations(InitialTension,X_o,**kwargs):
    """
	This function returns initial muscle lengths and muscle activations for a given pretensioning level, as derived from (***insert file_name here for scratchwork***) for the system that starts from rest. (Ex. pendulum_eqns.reference_trajectories._01).

	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	**kwargs

	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	1) Seed - Can see the random tension generated. When FixedInitialTension is provided, this seed will apply only to the initial conditions for activation and muscle length.

    2) PlotBool - Must be either True or False. Default is False. Will plot all possible initial muscle lengths and activations for a given pretensioning level.

    3) InitialTensionAcceleration - must be a numpy array of shape (2,). Default is set to the value generated from zero IC's. If using different reference trajectory, set InitialAngularAcceleration to d2r(0) (See below).

    4) InitialAngularAcceleration - must be either a numpy.float64, float, or int. Default is set to 0 to simulate starting from rest. Choice of reference trajectory *should* not matter as it is either 0 or d2r(0) (either by convention or by choice).

    5) InitialAngularSnap - must be either a numpy.float64, float, or int. Default is set to 0 to simulate starting from rest. Choice of reference trajectory *should* not matter as it is either 0 or d4r(0) (either by convention or by choice).

    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	"""
    PlotBool = kwargs.get("PlotBool",False)
    assert type(PlotBool)==bool,"PlotBool must be a boolean. Default is False."

    InitialAngularAcceleration = kwargs.get(
                "InitialAngularAcceleration",
                0
                ) # 0 or d2r(0)
    assert str(type(InitialAngularAcceleration)) in ["<class 'float'>","<class 'int'>","<class 'numpy.float64'>"], "InitialAngularAcceleration must be either a float or an int."

    InitialAngularSnap = kwargs.get(
                "InitialAngularSnap",
                0
                ) # 0 or d4r(0)
    assert str(type(InitialAngularSnap)) in ["<class 'float'>","<class 'int'>","<class 'numpy.float64'>"], "InitialAngularSnap must be either a float or an int."

    InitialTensionAcceleration = kwargs.get(
                "InitialTensionAcceleration",
                return_initial_tension_acceleration(
                    InitialTension,
                    X_o,
                    InitialAngularAcceleration=InitialAngularAcceleration,
                    InitialAngularSnap=InitialAngularSnap
                    )
                )
    assert np.shape(InitialTensionAcceleration)==(2,) \
    		and str(type(InitialTensionAcceleration))=="<class 'numpy.ndarray'>", \
    	"InitialTensionAcceleration must be a numpy array of shape (2,)"


    a_MTU1_o = np.sign(-r1(X_o[0]))*(
    	InitialAngularAcceleration
    	* np.sqrt(dr1_dθ(X_o[0])**2 + r1(X_o[0])**2)
    	+
    	X_o[1]**2
    	* dr1_dθ(X_o[0])
    	* (d2r1_dθ2(X_o[0]) + r1(X_o[0]))
    	/ np.sqrt(dr1_dθ(X_o[0])**2 + r1(X_o[0])**2)
    	)
    a_MTU2_o = np.sign(-r2(X_o[0]))*(
    	InitialAngularAcceleration
    	* np.sqrt(dr2_dθ(X_o[0])**2 + r2(X_o[0])**2)
    	+
    	X_o[1]**2
    	* dr2_dθ(X_o[0])
    	* (d2r2_dθ2(X_o[0]) + r2(X_o[0]))
    	/ np.sqrt(dr2_dθ(X_o[0])**2 + r2(X_o[0])**2)
    	)

    L1_UB = lo1*L_CE_max_1*(
    		k_1*np.log(
    				np.exp(
    					(m1*InitialTensionAcceleration[0]
    					+ (F_MAX1*cT/lTo1)
    						* (1-np.exp(-InitialTension[0]/(F_MAX1*cT*kT)))
    						* (c3*InitialTension[0]
    							- m1*a_MTU1_o
    						)
    					)
    					/ (F_MAX1*c3**2
    						*c_1*k_1
    						*(F_MAX1*cT/lTo1)
    						*(1-np.exp(-InitialTension[0]/(F_MAX1*cT*kT)))
    					)
    				)
    				- 1
    			)
    			+ Lr1
    		)
    L2_UB = lo2*L_CE_max_2*(
    		k_1*np.log(
    				np.exp(
    					(m2*InitialTensionAcceleration[1]
    					+ (F_MAX2*cT/lTo2)
    						* (1-np.exp(-InitialTension[1]/(F_MAX2*cT*kT)))
    						* (c4*InitialTension[1]
    							- m2*a_MTU2_o
    						)
    					)
    					/ (F_MAX2*c4**2
    						*c_1*k_1
    						*(F_MAX2*cT/lTo2)
    						*(1-np.exp(-InitialTension[1]/(F_MAX2*cT*kT)))
    					)
    				)
    				- 1
    			)
    			+ Lr1
    		)

    L1_LB = 0.5*lo1
    if L1_UB > 1.5*lo1:
    	L1_UB = 1.5*lo1
    L1 = np.linspace(L1_LB, L2_UB, 1001)
    # mu1, sigma1 = lo1, 0.1*lo1
    # L1 = np.array(list(sorted(np.random.normal(mu1, sigma1, 1001))))
    U1 = (m1*InitialTensionAcceleration[0]
    		+ (F_MAX1*cT/lTo1)
    			* (1-np.exp(-InitialTension[0]/(F_MAX1*cT*kT)))
    			* (c3*InitialTension[0]
    				- m1*a_MTU1_o
    				- F_MAX1*c3**3
    					*c_1*k_1
    					*np.log(np.exp((L1/(lo1*L_CE_max_1) - Lr1)/k_1)+1)
    				)
    	) \
    	/ (
    		F_MAX1*c3**2
    		*(F_MAX1*cT/lTo1)
    		*(1-np.exp(-InitialTension[0]/(F_MAX1*cT*kT)))
    		*np.exp(-(abs((L1-lo1)/(lo1*ω))**ρ))
    	)
    # U1 = (
    # 	InitialTension[0][0]/(F_MAX1*np.cos(α1))
    #     - c_1*k_1*np.log(np.exp((L1/(lo1*L_CE_max_1) - Lr1)/k_1)+1)
    #     ) / (np.exp(-(abs((L1-lo1)/(lo1*ω))**ρ)))

    L2_LB = 0.5*lo2
    if L2_UB > 1.5*lo2:
    	L2_UB = 1.5*lo2
    L2 = np.linspace(L2_LB, L2_UB, 1001)
    # mu2, sigma2 = lo2, 0.1*lo2
    # L2 = np.array(list(sorted(np.random.normal(mu2, sigma2, 1001))))
    U2 = (m2*InitialTensionAcceleration[1]
    		+ (F_MAX2*cT/lTo2)
    			* (1-np.exp(-InitialTension[1]/(F_MAX2*cT*kT)))
    			* (c4*InitialTension[1]
    				- m2*a_MTU2_o
    				- F_MAX2*c4**3
    					*c_1*k_1
    					*np.log(np.exp((L2/(lo2*L_CE_max_2) - Lr1)/k_1)+1)
    				)
    	) \
    	/ (
    		F_MAX2*c4**2
    		*(F_MAX2*cT/lTo2)
    		*(1-np.exp(-InitialTension[1]/(F_MAX2*cT*kT)))
    		*np.exp(-(abs((L2-lo2)/(lo2*ω))**ρ))
    	)
    # U2 = (
    # 	InitialTension[1][0]/(F_MAX2*np.cos(α2))
    #     - c_1*k_1*np.log(np.exp((L2/(lo2*L_CE_max_2) - Lr1)/k_1)+1)
    #     ) / (np.exp(-(abs((L2-lo2)/(lo2*ω))**ρ)))

    if PlotBool == True:
    	plt.figure(figsize=(10,8))
    	plt.title(r"Viable Initial $l_{m,1}$ and $u_{1}$ Values")
    	plt.xlabel(r"$l_{m,1}$ (m)",fontsize=14)
    	plt.ylabel(r"$u_{1}$",fontsize=14)
    	plt.scatter(L1,U1)
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
    	plt.scatter(L2,U2)
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
    return(L1,U1,L2,U2)
#
# def return_random_initial_muscle_lengths_and_activations(InitialTension,**kwargs):
#     """
# 	This function returns initial conditions for the system that starts from rest. (Ex. pendulum_eqns.reference_trajectories._01).
#
# 	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# 	**kwargs
#
# 	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# 	1) Seed - Can see the random tension generated. When FixedInitialTension is provided, this seed will apply only to the initial conditions for activation and muscle length.
#
#     2) PlotBool - Must be either True or False. Default is False. Will plot all possible initial muscle lengths and activations for a given pretensioning level.
#
#     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 	"""
#
#     Seed = kwargs.get("Seed",None)
# 	assert type(Seed) in [float,int] or Seed is None, "Seed must be a float or an int or None."
# 	np.random.seed(Seed)
#
# 	PlotBool = kwargs.get("PlotBool",False)
# 	assert type(PlotBool)==bool,"PlotBool must be a boolean. Default is False."
#
# 	# L1 = np.linspace(0.5*lo1, 1.5*lo1, 1001)
# 	mu1, sigma1 = lo1, 0.1*lo1
# 	L1 = np.array(list(sorted(np.random.normal(mu1, sigma1, 1001))))
# 	U1 = (
# 		InitialTension[0][0]/(F_MAX1*np.cos(α1))
# 	    - c_1*k_1*np.log(np.exp((L1/(lo1*L_CE_max_1) - Lr1)/k_1)+1)
# 	    ) / (np.exp(-(abs((L1-lo1)/(lo1*ω))**ρ)))
#
# 	# L2 = np.linspace(0.5*lo2, 1.5*lo2, 1001)
# 	mu2, sigma2 = lo2, 0.1*lo2
# 	L2 = np.array(list(sorted(np.random.normal(mu2, sigma2, 1001))))
# 	U2 = (
# 		InitialTension[1][0]/(F_MAX2*np.cos(α2))
# 	    - c_1*k_1*np.log(np.exp((L2/(lo2*L_CE_max_2) - Lr1)/k_1)+1)
# 	    ) / (np.exp(-(abs((L2-lo2)/(lo2*ω))**ρ)))
#
# 	if PlotBool == True:
# 		plt.figure(figsize=(10,8))
# 		plt.title(r"Viable Initial $l_{m,1}$ and $u_{1}$ Values")
# 		plt.xlabel(r"$l_{m,1}$ (m)",fontsize=14)
# 		plt.ylabel(r"$u_{1}$",fontsize=14)
# 		plt.scatter(L1,U1)
# 		plt.plot([lo1,lo1],[0,1],'0.70',linestyle='--')
# 		plt.gca().set_ylim((0,1))
# 		plt.gca().set_xticks(
# 			[0.25*lo1,
# 			0.5*lo1,
# 			0.75*lo1,
# 			lo1,
# 			1.25*lo1,
# 			1.5*lo1,
# 			1.75*lo1]
# 			)
# 		plt.gca().set_xticklabels(
# 			["",
# 			r"$\frac{1}{2}$ $l_{o,2}$",
# 			"",
# 			r"$l_{o,2}$",
# 			"",
# 			r"$\frac{3}{2}$ $l_{o,2}$",
# 			""],
# 			fontsize=12)
#
# 		plt.figure(figsize=(10,8))
# 		plt.title(r"Viable Initial $l_{m,2}$ and $u_{2}$ Values")
# 		plt.xlabel(r"$l_{m,2}$ (m)",fontsize=14)
# 		plt.ylabel(r"$u_{2}$",fontsize=14)
# 		plt.scatter(L2,U2)
# 		plt.plot([lo2,lo2],[0,1],'0.70',linestyle='--')
# 		plt.gca().set_ylim((0,1))
# 		plt.gca().set_xticks(
# 			[0.25*lo2,
# 			0.5*lo2,
# 			0.75*lo2,
# 			lo2,
# 			1.25*lo2,
# 			1.5*lo2,
# 			1.75*lo2]
# 			)
# 		plt.gca().set_xticklabels(
# 			["",
# 			r"$\frac{1}{2}$ $l_{o,2}$",
# 			"",
# 			r"$l_{o,2}$",
# 			"",
# 			r"$\frac{3}{2}$ $l_{o,2}$",
# 			""],
# 			fontsize=12)
#
# 		plt.show()
# 	return(L1,U1,L2,U2)

def find_viable_initial_values(**kwargs):
	"""
	This function returns initial conditions for the system that starts from rest. (Ex. pendulum_eqns.reference_trajectories._01)

	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	**kwargs (Parent)

	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	1) FixedInitialTension - Must be a (2,1) numpy.ndarray. Run find_initial_tension outside of the loop for a given seed and then feed it through the pipeline.

	2) ReturnAll - Can return all initial values for a given tension level. Will be fed through to return_random_initial_muscle_lengths_and_activations.

	3) Seed - Can see the random tension generated. When FixedInitialTension is provided, this seed will apply only to the initial conditions for activation and muscle length.

    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	**kwargs (Passed to return_initial_tension(X_o,**kwargs))

    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	4) InitialAngularAcceleration - must be a float or an int. Default is 0 (starting from rest).

	5) Seed - must be a float or an int. Default is None (seeded by current time).

    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	**kwargs (Passed to return_random_initial_muscle_lengths_and_activations(T_o,X_o,**kwargs))

	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	6) Seed - Can see the random tension generated. When FixedInitialTension is provided, this seed will apply only to the initial conditions for activation and muscle length.

    7) PlotBool - Must be either True or False. Default is False. Will plot all possible initial muscle lengths and activations for a given pretensioning level.

    8) InitialTensionAcceleration - must be a numpy array of shape (2,). Default is set to the value generated from joint angle IC's for pendulum_eqns.reference_trajectories._01.py (B+A,0,-Aw²,0,Aw⁴,etc.). If using a different reference trajectory, it would be best to either set all joint angle IC's to zero (like for pendulum_eqns.reference_trajectories._02.py) or to derive the value and pass it through the kwargs here (NOTE: must also consider how changing angular IC's effect other state derivative IC's).

    9) InitialAngularAcceleration - must be either a numpy.float64, float, or int. Default is set to 0 to simulate starting from rest. Choice of reference trajectory *should* not matter as it is either 0 or d2r(0) (either by convention or by choice).

    10) InitialAngularSnap - must be either a numpy.float64, float, or int. Default is set to 0 to simulate starting from rest. Choice of reference trajectory *should* not matter as it is either 0 or d4r(0) (either by convention or by choice).

    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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

	X_o = np.array([r(0),dr(0)])
	if FixedInitialTension is None:
		T = return_initial_tension(X_o,**kwargs)
	else:
		T = FixedInitialTension
	L1,U1,L2,U2 = return_random_initial_muscle_lengths_and_activations(T,X_o,**kwargs)
	rand_index = np.random.choice(len(L1),2)

	if ReturnAll == False:
		return(
			T,
			np.array([L1[rand_index[0]],L2[rand_index[1]]]),
			np.array([U1[rand_index[0]],U2[rand_index[1]]])
			)
	else:
		return(T,L1,L2,U1,U2)
