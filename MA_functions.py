from muscle_settings import *

def MA_function(Parameters,θ_PS=None):
	"""
	Note:

	Angles should be a number if Coefficients has a length of 5, or a list of length 2 when the Coefficients have lengths 16 or 18. Angles[0] will be the PRIMARY ANGLE for the DOF being considered while Angles[1] will be the secondary angle.

	Notes:

	threshold is only needed for Pigeon or Ramsay; 2009 MA functions that are invalid outside of a given value. Must be either None (default) or the radian value of the threshold.

	eq is only needed for Ramsay; 2009 (Pigeon has one quintic polynomial). eq must be either 1, 2, or 3, with list length requirements of 5, 16, or 18, respectively.
	"""

	Parameters = return_primary_source(Parameters)
	assert str(type(Parameters))=="<class 'muscle_settings.MA_Settings'>", "Parameters are not in correct namedtuple form. Should be <class 'muscle_settings.MA_Settings'> instead of " + str(type(Parameters))
	if θ_PS is None:
		θ_PS = np.pi
	else:
		assert type(θ_PS)==float, "θ_PS must be a float."
	src = Parameters.Source
	Coefficients = unit_conversion(Parameters)
	eq = Parameters.Equation_Number
	threshold = Parameters.Threshold

	assert type(src) == str, "src must be a str."
	assert src.capitalize() in ['Ramsay; 2009','Pigeon; 1996','Kuechle; 1997','Holzbaur; 2005', 'Est'], "src must be either Ramsay; 2009, Pigeon or Est (Estimate)."

	'''
	Note:
	For Kuechle and Holzbaur, where estimates or average MA were given, the format should be [MA,0,0,0,0,0] such that the function returns a constant MA function (See matrix multiplication below).
	'''

	if src.capitalize() in ['Pigeon; 1996', 'Kuechle; 1997', 'Holzbaur; 2005']:
		assert len(Coefficients)==6, 'For Pigeon (1996) the list of Coefficients must be 6 elements long. Insert zeros (0) for any additional empty coefficients.'
		MomentArm = lambda θ: (np.matrix(Coefficients,dtype='float64')\
										*np.matrix([1,θ,θ**2,θ**3,θ**4,θ**5]).T)[0,0]
	elif src.capitalize() == 'Est':
		MomentArm = lambda θ: np.array(Coefficients,dtype='float64')
	else: #src.capitalize() == 'Ramsay; 2009'
		assert type(Coefficients) == list, "Coefficients must be a list."
		assert len(Coefficients) in [5,16,18], "Coefficients as a list must be of length 5, 16, or 18."
		assert eq in [1,2,3], "eq must be either 1, 2, or 3 when using Ramsay; 2009 (2009)."
		if eq == 1:
			assert len(Coefficients) == 5, "For Eq. 1, Coefficients must be 5 elements long."
			MomentArm = lambda θ: \
					(np.matrix(Coefficients,dtype='float64')*\
						np.matrix([1,θ,θ**2,θ**3,θ**4]).T)[0,0]
		elif eq == 2:
			assert len(Coefficients)==16, "For Eq. 2, Coefficients must be 16 elements long."
			MomentArm = lambda θ: \
					(np.matrix(Coefficients,dtype='float64')*\
						np.matrix([1, θ, θ_PS, θ*θ_PS, θ**2, \
									θ_PS**2, (θ**2)*θ_PS, θ*(θ_PS**2), \
									(θ**2)*(θ_PS**2), θ**3, θ_PS**3, \
									(θ**3)*θ_PS, θ*(θ_PS**3), \
									(θ**3)*(θ_PS**2), (θ**2)*(θ_PS**3), \
									(θ**3)*(θ_PS**3)]).T)[0, 0]
		else: # eq == 3
			assert len(Coefficients)==18, "For Eq. 3, Coefficients must be 18 elements long."
			MomentArm = lambda θ: \
					(np.matrix(Coefficients,dtype='float64')*\
						np.matrix([1, θ, θ_PS, θ*θ_PS, θ**2, \
									θ_PS**2, (θ**2)*θ_PS, θ*(θ_PS**2), (θ**2)*(θ_PS**2), \
									θ**3, (θ**3)*θ_PS, (θ**3)*(θ_PS**2), \
									θ**4, (θ**4)*θ_PS, (θ**4)*(θ_PS**2),  \
									θ**5, (θ**5)*θ_PS, (θ**5)*(θ_PS**2)]).T)[0, 0]
	if threshold is None:
		return(MomentArm)
	else:
		assert type(threshold) in [int,float], "threshold must be a number."
		PiecewiseMomentArm = lambda θ:\
					np.piecewise(θ,[θ<threshold,θ>=threshold],\
									[MomentArm(θ),MomentArm(threshold)])
		return(PiecewiseMomentArm)
def MA_deriv(Parameters,θ_PS=None):
	"""
	Note:

	Angles should be a number if Coefficients has a length of 5, or a list of length 2 when the Coefficients have lengths 16 or 18. Angles[0] will be the PRIMARY ANGLE for the DOF being considered while Angles[1] will be the secondary angle.

	Notes:

	threshold is only needed for Pigeon or Ramsay; 2009 MA functions that are invalid outside of a given value. Must be either None (default) or the radian value of the threshold.

	eq is only needed for Ramsay; 2009 (Pigeon has one quintic polynomial). eq must be either 1, 2, or 3, with list length requirements of 5, 16, or 18, respectively.
	"""

	Parameters = return_primary_source(Parameters)
	assert str(type(Parameters))=="<class 'muscle_settings.MA_Settings'>", "Parameters are not in correct namedtuple form."
	if θ_PS is None:
		θ_PS = np.pi
	else:
		assert type(θ_PS)==float, "θ_PS must be a float."
	src = Parameters.Source
	Coefficients = unit_conversion(Parameters)
	eq = Parameters.Equation_Number
	threshold = Parameters.Threshold

	assert type(src) == str, "src must be a str."
	assert src.capitalize() in ['Ramsay; 2009','Pigeon; 1996','Kuechle; 1997','Holzbaur; 2005', 'Est'], "src must be either Ramsay; 2009, Pigeon or Est (Estimate)."

	'''
	Note:
	For Kuechle and Holzbaur, where estimates or average MA were given, the format should be [MA,0,0,0,0,0] such that the function returns a constant MA function (See matrix multiplication below).
	'''

	if src.capitalize() in ['Pigeon; 1996', 'Kuechle; 1997', 'Holzbaur; 2005']:
		assert len(Coefficients)==6, 'For Pigeon (1996) the list of Coefficients must be 6 elements long. Insert zeros (0) for any additional empty coefficients.'
		"""
		(d/dθ) [MomentArm] = (np.matrix(Coefficients,dtype='float64')\
						*np.matrix([(d/dθ)[1],(d/dθ)[θ],(d/dθ)[θ**2],\
										(d/dθ)[θ**3],(d/dθ)[θ**4],(d/dθ)[θ**5]]).T)[0,0]
		"""
		Derivative = lambda θ: (np.matrix(Coefficients,dtype='float64')\
						*np.matrix([0,1,(2*θ),(3*θ**2),(4*θ**3),(5*θ**4)]).T)[0,0]
	elif src.capitalize() == 'Est':
		"""
		(d/dθ)[MomentArm] = np.array((d/dθ)[Coefficients],dtype='float64')
		"""
		Derivative = lambda θ:  0
	else: #src.capitalize() == 'Ramsay; 2009'
		assert type(Coefficients) == list, "Coefficients must be a list."
		assert len(Coefficients) in [5,16,18], "Coefficients as a list must be of length 5, 16, or 18."
		assert eq in [1,2,3], "eq must be either 1, 2, or 3 when using Ramsay; 2009 (2009)."
		if eq == 1:
			assert len(Coefficients) == 5, "For Eq. 1, Coefficients must be 5 elements long."
			"""
			(d/dθ)[MomentArm] = (np.matrix(Coefficients,dtype='float64')\
									*np.matrix([	(d/dθ)[1],		(d/dθ)[θ],\
													(d/dθ)[θ**2],	(d/dθ)[θ**3],\
													(d/dθ)[θ**4]]).T)[0,0]
			"""
			Derivative = lambda θ: (np.matrix(Coefficients,dtype='float64')\
							*np.matrix([	0,				1,\
											(2*θ),			(3*θ**2),\
											(4*θ**3)]).T)[0,0]
		elif eq == 2:
			"""
			This is only good for this ReferenceTracking Ex where PS is fixed. Derivative is only wrt one DOF.
			"""
			assert len(Coefficients)==16, "For Eq. 2, Coefficients must be 16 elements long."
			"""
			(d/dθ)[MomentArm] = \
					(np.matrix(Coefficients,dtype='float64')*\
							np.matrix([(d/dθ)[1], 					(d/dθ)[θ], \
										(d/dθ)[θ_PS], 				(d/dθ)[θ*θ_PS],\
										(d/dθ)[θ**2],				(d/dθ)[θ_PS**2],\
										(d/dθ)[(θ**2)*θ_PS],		(d/dθ)[θ*(θ_PS**2)],\
										(d/dθ)[(θ**2)*(θ_PS**2)],	(d/dθ)[θ**3], \
										(d/dθ)[θ_PS**3], 			(d/dθ)[(θ**3)*θ_PS],\
										(d/dθ)[θ*(θ_PS**3)],		(d/dθ)[(θ**3)*(θ_PS**2)],\
										(d/dθ)[(θ**2)*(θ_PS**3)],	(d/dθ)[(θ**3)*(θ_PS**3)]\
										]).T)[0, 0]
			"""
			Derivative = lambda θ: (np.matrix(Coefficients,dtype='float64')*\
							np.matrix([ 0, 					1,\
										0,  				θ_PS, \
										(2*θ), 				0, \
										(2*θ)*θ_PS, 		(θ_PS**2), 	\
										(2*θ)*(θ_PS**2), 	(3*θ**2), \
										0, 					(3*θ**2)*θ_PS, \
										(θ_PS**3), 			(3*θ**2)*(θ_PS**2), \
										(2*θ)*(θ_PS**3),	(3*θ**2)*(θ_PS**3)]).T)[0, 0]
		else: # eq == 3
			assert len(Coefficients)==18, "For Eq. 3, Coefficients must be 18 elements long."
			"""
			(d/dθ)[MomentArm] = \
					(np.matrix(Coefficients,dtype='float64')*\
							np.matrix([(d/dθ)[1], 					(d/dθ)[θ], \
										(d/dθ)[θ_PS], 				(d/dθ)[θ*θ_PS],\
										(d/dθ)[θ**2], 				(d/dθ)[θ_PS**2],\
										(d/dθ)[(θ**2)*θ_PS], 		(d/dθ)[θ*(θ_PS**2)],\
										(d/dθ)[(θ**2)*(θ_PS**2)], 	(d/dθ)[θ**3],\
										(d/dθ)[(θ**3)*θ_PS], 		(d/dθ)[(θ**3)*(θ_PS**2)],\
										(d/dθ)[θ**4], 				(d/dθ)[(θ**4)*θ_PS],\
										(d/dθ)[(θ**4)*(θ_PS**2)],  	(d/dθ)[θ**5],\
										(d/dθ)[(θ**5)*θ_PS], 		(d/dθ)[(θ**5)*(θ_PS**2)\
										]]).T)[0, 0]
			"""
			Derivative = lambda θ: (np.matrix(Coefficients,dtype='float64')*\
								np.matrix([	0, 					1,\
											0, 					θ_PS,\
										  	(2*θ),				0,\
									 		(2*θ)*θ_PS, 		(θ_PS**2),\
									  		(2*θ)*(θ_PS**2),	(3*θ**2),\
									 		(3*θ**2)*θ_PS, 		(3*θ**2)*(θ_PS**2),\
											(4*θ**3), 			(4*θ**3)*θ_PS,\
									 		(4*θ**3)*(θ_PS**2),	(5*θ**4),\
									 		(5*θ**4)*θ_PS, 		(5*θ**4)*(θ_PS**2)]).T)[0, 0]
	if threshold is None:
		return(Derivative)
	else:
		assert type(threshold) in [int,float], "threshold must be a number."
		PiecewiseDerivative = lambda θ:\
									np.piecewise(θ,[θ<threshold,θ>=threshold],\
													[Derivative(θ),0])
		return(PiecewiseDerivative)
def MA_2nd_deriv(Parameters,θ_PS=None):
	"""
	Note:

	Angles should be a number if Coefficients has a length of 5, or a list of length 2 when the Coefficients have lengths 16 or 18. Angles[0] will be the PRIMARY ANGLE for the DOF being considered while Angles[1] will be the secondary angle.

	Notes:

	threshold is only needed for Pigeon or Ramsay; 2009 MA functions that are invalid outside of a given value. Must be either None (default) or the radian value of the threshold.

	eq is only needed for Ramsay; 2009 (Pigeon has one quintic polynomial). eq must be either 1, 2, or 3, with list length requirements of 5, 16, or 18, respectively.
	"""

	import numpy as np
	Parameters = return_primary_source(Parameters)
	assert str(type(Parameters))=="<class 'muscle_settings.MA_Settings'>", "Parameters are not in correct namedtuple form."
	if θ_PS is None:
		θ_PS = np.pi
	else:
		assert type(θ_PS)==float, "θ_PS must be a float."
	src = Parameters.Source
	Coefficients = unit_conversion(Parameters)
	eq = Parameters.Equation_Number
	threshold = Parameters.Threshold

	assert type(src) == str, "src must be a str."
	assert src.capitalize() in ['Ramsay; 2009','Pigeon; 1996','Kuechle; 1997','Holzbaur; 2005', 'Est'], "src must be either Ramsay; 2009, Pigeon or Est (Estimate)."

	'''
	Note:
	For Kuechle and Holzbaur, where estimates or average MA were given, the format should be [MA,0,0,0,0,0] such that the function returns a constant MA function (See matrix multiplication below).
	'''

	if src.capitalize() in ['Pigeon; 1996', 'Kuechle; 1997', 'Holzbaur; 2005']:
		assert len(Coefficients)==6, 'For Pigeon (1996) the list of Coefficients must be 6 elements long. Insert zeros (0) for any additional empty coefficients.'
		"""
		(d²/dθ²) [MomentArm] \
			= (np.matrix(Coefficients,dtype='float64')\
						*np.matrix([(d²/dθ²)[1],		(d²/dθ²)[θ],\
									(d²/dθ²)[θ**2],		(d²/dθ²)[θ**3],\
									(d²/dθ²)[θ**4],		(d²/dθ²)[θ**5]]).T)[0,0]

			= (d/dθ)[Derivative]
			= (np.matrix(Coefficients,dtype='float64')\
						*np.matrix([(d/dθ)[0],			(d/dθ)[1],\
									(d/dθ)[(2*θ)],		(d/dθ)[(3*θ**2)],\
									(d/dθ)[(4*θ**3)],	(d/dθ)[(5*θ**4)]]).T)[0,0]
		"""

		SecondDerivative = lambda θ: (np.matrix(Coefficients,dtype='float64')\
											*np.matrix([0,			0,\
														2,			6*θ,\
														(12*θ**2),	(20*θ**3)]).T)[0,0]
	elif src.capitalize() == 'Est':
		"""
		(d²/dθ²)[MomentArm] = np.array((d²/dθ²)[Coefficients],dtype='float64')
		"""
		# (d/dθ)[Derivative] = lambda θ:  (d/dθ)[0]
		SecondDerivative = lambda θ:  0
	else: #src.capitalize() == 'Ramsay; 2009'
		assert type(Coefficients) == list, "Coefficients must be a list."
		assert len(Coefficients) in [5,16,18], "Coefficients as a list must be of length 5, 16, or 18."
		assert eq in [1,2,3], "eq must be either 1, 2, or 3 when using Ramsay; 2009 (2009)."
		if eq == 1:
			assert len(Coefficients) == 5, "For Eq. 1, Coefficients must be 5 elements long."
			"""
			(d²/dθ²)[MomentArm] \
				= (np.matrix(Coefficients,dtype='float64')\
							*np.matrix([(d²/dθ²)[1],		(d²/dθ²)[θ],\
										(d²/dθ²)[θ**2],		(d²/dθ²)[θ**3],\
										(d²/dθ²)[θ**4]]).T)[0,0]

				= (d/dθ)[Derivative] \
				= (np.matrix(Coefficients,dtype='float64')\
							*np.matrix([(d/dθ)[0],			(d/dθ)[1],\
										(d/dθ)[(2*θ)],  	(d/dθ)[(3*θ**2)],\
										(d/dθ)[(4*θ**3)]]).T)[0,0]
			"""
			SecondDerivative = lambda θ: \
						(np.matrix(Coefficients,dtype='float64')\
								*np.matrix([0,			0,\
											2, 			6*θ,\
											(12*θ**2)				]).T)[0,0]
		elif eq == 2:
			"""
			This is only good for this ReferenceTracking Ex where PS is fixed. Derivative is only wrt one DOF.
			"""
			assert len(Coefficients)==16, "For Eq. 2, Coefficients must be 16 elements long."
			"""
			(d²/dθ²)[MomentArm] \
				= (np.matrix(Coefficients,dtype='float64')*\
						np.matrix([	(d²/dθ²)[1], 				(d²/dθ²)[θ], \
									(d²/dθ²)[θ_PS], 			(d²/dθ²)[θ*θ_PS],\
									(d²/dθ²)[θ**2],				(d²/dθ²)[θ_PS**2],\
									(d²/dθ²)[(θ**2)*θ_PS],		(d²/dθ²)[θ*(θ_PS**2)],\
									(d²/dθ²)[(θ**2)*(θ_PS**2)],	(d²/dθ²)[θ**3], \
									(d²/dθ²)[θ_PS**3], 			(d²/dθ²)[(θ**3)*θ_PS],\
									(d²/dθ²)[θ*(θ_PS**3)],		(d²/dθ²)[(θ**3)*(θ_PS**2)],\
									(d²/dθ²)[(θ**2)*(θ_PS**3)],	(d²/dθ²)[(θ**3)*(θ_PS**3)]\
									]).T)[0, 0]

				= (d/dθ)[Derivative] \
				= (np.matrix(Coefficients,dtype='float64')*\
						np.matrix([ (d/dθ)[0], 					(d/dθ)[1], \
									(d/dθ)[0], 					(d/dθ)[θ_PS], \
									(d/dθ)[(2*θ)], 				(d/dθ)[0], \
									(d/dθ)[(2*θ)*θ_PS], 		(d/dθ)[(θ_PS**2)], \
									(d/dθ)[(2*θ)*(θ_PS**2)], 	(d/dθ)[(3*θ**2)], \
									(d/dθ)[0], 					(d/dθ)[(3*θ**2)*θ_PS], \
									(d/dθ)[(θ_PS**3)], 			(d/dθ)[(3*θ**2)*(θ_PS**2)], \
									(d/dθ)[(2*θ)*(θ_PS**3)],	(d/dθ)[(3*θ**2)*(θ_PS**3)]\
									]).T)[0, 0]
			"""
			SecondDerivative = lambda θ: \
					(np.matrix(Coefficients,dtype='float64')*\
						np.matrix([ 0, 					0,\
						 			0,					0,\
									2,					0,\
				 					2*θ_PS,				0,\
						 			2*(θ_PS**2),		(6*θ),\
									0,					(6*θ)*θ_PS,\
									0,					(6*θ)*(θ_PS**2),\
									2*(θ_PS**3),		(6*θ)*(θ_PS**3)\
									]).T)[0, 0]
		else: # eq == 3
			assert len(Coefficients)==18, "For Eq. 3, Coefficients must be 18 elements long."
			"""
			(d²/dθ²)[MomentArm] \
				= (np.matrix(Coefficients,dtype='float64')*\
						np.matrix([(d²/dθ²)[1], 					(d²/dθ²)[θ], \
									(d²/dθ²)[θ_PS], 				(d²/dθ²)[θ*θ_PS],\
									(d²/dθ²)[θ**2], 				(d²/dθ²)[θ_PS**2],\
									(d²/dθ²)[(θ**2)*θ_PS], 			(d²/dθ²)[θ*(θ_PS**2)],\
									(d²/dθ²)[(θ**2)*(θ_PS**2)], 	(d²/dθ²)[θ**3],\
									(d²/dθ²)[(θ**3)*θ_PS], 			(d²/dθ²)[(θ**3)*(θ_PS**2)],\
									(d²/dθ²)[θ**4], 				(d²/dθ²)[(θ**4)*θ_PS],\
									(d²/dθ²)[(θ**4)*(θ_PS**2)], 	(d²/dθ²)[θ**5],\
									(d²/dθ²)[(θ**5)*θ_PS], 			(d²/dθ²)[(θ**5)*(θ_PS**2)\
									]]).T)[0, 0]

				= (d/dθ)[Derivative] \
				= (np.matrix(Coefficients,dtype='float64')*\
						np.matrix([	(d/dθ)[0], 					(d/dθ)[1],\
									(d/dθ)[0], 					(d/dθ)[θ_PS],\
								  	(d/dθ)[(2*θ)],				(d/dθ)[0],\
							 		(d/dθ)[(2*θ)*θ_PS], 		(d/dθ)[(θ_PS**2)],\
							  		(d/dθ)[(2*θ)*(θ_PS**2)],	(d/dθ)[(3*θ**2)],\
							 		(d/dθ)[(3*θ**2)*θ_PS], 		(d/dθ)[(3*θ**2)*(θ_PS**2)],\
									(d/dθ)[(4*θ**3)], 			(d/dθ)[(4*θ**3)*θ_PS],\
							 		(d/dθ)[(4*θ**3)*(θ_PS**2)],	(d/dθ)[(5*θ**4)],\
							 		(d/dθ)[(5*θ**4)*θ_PS], 		(d/dθ)[(5*θ**4)*(θ_PS**2)]\
									]).T)[0, 0]
			"""
			SecondDerivative = lambda θ: \
					(np.matrix(Coefficients,dtype='float64')*\
							np.matrix([	0, 						0,\
										0, 						0,\
										2,						0,\
										2*θ_PS, 				0,\
										2*(θ_PS**2),			(6*θ),\
										(6*θ)*θ_PS, 			(6*θ)*(θ_PS**2),\
										(12*θ**2), 				(12*θ**2)*θ_PS,\
										(12*θ**2)*(θ_PS**2),	(20*θ**3),\
										(20*θ**3)*θ_PS, 		(20*θ**3)*(θ_PS**2)\
										]).T)[0, 0]
	if threshold is None:
		return(SecondDerivative)
	else:
		assert type(threshold) in [int,float], "threshold must be a number."
		PiecewiseSecondDerivative = lambda θ:\
									np.piecewise(θ,[θ<threshold,θ>=threshold],\
													[SecondDerivative(θ),0])
		return(PiecewiseSecondDerivative)
def return_MA_matrix_functions(
		AllMuscleSettings,ReturnMatrixFunction=False,θ_PS=None):
	"""returns an (n,m) matrix when n is the number of muscles and m is the number of DOFS. We chose to return R.T because this is commonly utilized in muscle velocity calculations.
	"""
	import numpy as np
	import sympy as sp

	MuscleList = AllMuscleSettings.keys()
	if ReturnMatrixFunction == False:
		R_Transpose = np.matrix([\
			[MA_function(AllMuscleSettings[muscle]["Shoulder MA"],θ_PS=θ_PS), \
			MA_function(AllMuscleSettings[muscle]["Elbow MA"],θ_PS=θ_PS)]\
				for muscle in MuscleList])
		dR_Transpose = np.matrix([\
			[MA_deriv(AllMuscleSettings[muscle]["Shoulder MA"],θ_PS=θ_PS), \
			MA_deriv(AllMuscleSettings[muscle]["Elbow MA"],θ_PS=θ_PS)]\
				for muscle in MuscleList])
		d2R_Transpose = np.matrix([\
			[MA_2nd_deriv(AllMuscleSettings[muscle]["Shoulder MA"],θ_PS=θ_PS), \
			MA_2nd_deriv(AllMuscleSettings[muscle]["Elbow MA"],θ_PS=θ_PS)]\
					for muscle in MuscleList])
	else:
		R_Transpose = lambda θ_SFE, θ_EFE: \
			np.matrix([[MA_function(AllMuscleSettings[muscle]["Shoulder MA"],θ_PS=θ_PS)(θ_SFE),\
						MA_function(AllMuscleSettings[muscle]["Elbow MA"],θ_PS=θ_PS)(θ_EFE)] \
							for muscle in MuscleList])
		dR_Transpose = lambda θ_SFE, θ_EFE, θ_PS: \
			np.matrix([[MA_deriv(AllMuscleSettings[muscle]["Shoulder MA"],θ_PS=θ_PS)(θ_SFE),\
						MA_deriv(AllMuscleSettings[muscle]["Elbow MA"],θ_PS=θ_PS)(θ_EFE)] \
							for muscle in MuscleList])
		d2R_Transpose = lambda θ_SFE, θ_EFE, θ_PS: \
			np.matrix([[MA_2nd_deriv(AllMuscleSettings[muscle]["Shoulder MA"],θ_PS=θ_PS)(θ_SFE),\
						MA_2nd_deriv(AllMuscleSettings[muscle]["Elbow MA"],θ_PS=θ_PS)(θ_EFE)] \
							for muscle in MuscleList])
	return(R_Transpose,dR_Transpose,d2R_Transpose)
