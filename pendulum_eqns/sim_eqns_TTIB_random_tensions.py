import numpy as np
import cvxopt
import matplotlib.pyplot as plt
from termcolor import cprint,colored
from danpy.sb import dsb,get_terminal_width
from pendulum_eqns.init_tendon_tension_controlled_model import *

N_seconds = 4
N = N_seconds*100 + 1
Time = np.linspace(0,N_seconds,N)
dt = Time[1]-Time[0]

def split_time_array(N,Time):
    correctNumberOfBins = False
    while correctNumberOfBins == False:
        timeArrayLengths = [int((len(Time)/N)*el) for el in np.random.normal(1,0.1,size=(N,))]
        if sum(timeArrayLengths)==len(Time):
            correctNumberOfBins=True
    indices = list(np.cumsum(timeArrayLengths))
    indices.append(0)
    indices = list(sorted(indices))

    splitTimeArrays = []
    for i in range(len(indices)-1):
        splitTimeArrays.append(Time[indices[i]:indices[i+1]])

    return(splitTimeArrays)

class Spline:
    """
    Initiate a class variable spline that has one break at x = x_break starting at x_initial and has
    the equation y = a + b*(x-x_o) + c*(x-x_o)**2 + d*(x-x_o)**3.

    pp_func()
    ~~~~~~~~~~~~~~~~~~~

    Takes in X array and outputs the piecewise polynomial associated with this spline.

    pp_deriv()
    ~~~~~~~~~~~~~~~~~~~

    Takes in X array and outputs the piecewise polynomial associated with the spline's derivative.

    pp_2deriv()
    ~~~~~~~~~~~~~~~~~~~

    Takes in X array and outputs the piecewise polynomial associated with the spline's second derivative.

    find_max_and_min()
    ~~~~~~~~~~~~~~~~~~~

    Takes in the min and max values for both x and y and will find the maximum y values of the piecewise
    polynomial. To do this, first we find the extrema point (find_extrema) by inputing the x values that
    set the derivate of the piecewise polynomial equal to zero (quadratic formula). Next we ensure that
    the zero values are in fact real (is_real). We then filter out the zeros that are not in the
    appropriate domains (is_in_appropriate_domain). To see if these values are maximum or minimum, we
    plug them back into the second derivative of the appropriate piecewise polynomial (second_deriv_is_neg()
    and second_deriv_is_pos(), respectively). Finally we determine the y value of these extrema by using
    the class function self.pp_func().

    is_initial_slope_positive()
    ~~~~~~~~~~~~~~~~~~~

    This takes in X and will check to make sure that for the first 2500 entries in X, that the derivative
    of the piecewise polynomial (pp_deriv()) will be positive. Make sure that X is at least 2500 in length.

    is_within_bounds()
    ~~~~~~~~~~~~~~~~~~~

    This checks to see if the maximum maximum value and the minimum mininum value calculated above will fall between y_min and y_max. This makes use of self.find_max_and_min()

    print_func()
    ~~~~~~~~~~~~~~~~~~~

    This function uses pprint() to return a printout of the piecewise polynomial f(x).

    return_parameterized_X()
    ~~~~~~~~~~~~~~~~~~~

    This function will return the parameterized x(t) and y(t) that follows path S along the curve y=f(x) subject to some tangential velocity profile dS/dt. This utilizes scipy.integrate.odeint to solve the time derivative to the path length function S = ∫(dx_dt)√(1 + f'(x)²)dt.

    return_parameterized_dX()
    ~~~~~~~~~~~~~~~~~~~

    This function takes in x(t) and returns dx/dt and dy/dt derived from the relationship of arc length S, its time derivative dS/dt, and the path f(x(t)) and its derivative with respect to x.

    return_parameterized_d2X()
    ~~~~~~~~~~~~~~~~~~~

    This function takes in x(t) and dx/dt and returns d²x/dt² and d²y/dt² derived from the relationship of arc length S, its first and second time derivatives, dS/dt and d²S/dt², respectively, and the path f(x(t)) and its first and second derivatives with respect to x.

    find_path_length()
    ~~~~~~~~~~~~~~~~~~~

    Calculates the path length of the curve y=f(x) from S = ∫(dx_dt)√(1 + f'(x)²)dt. This is needed in order to describe the minimum jerk criterion tangential velocity equation, dS/dt.

    dS_dt()
    ~~~~~~~~~~~~~~~~~~~

    Returns the minimum jerk criterion tangential velocity equation, dS/dt, given by:

    						dS/dt = S*(30t² - 60t³ + 30t⁴)

    Where S is found from S = self.find_path_length()

    d2S_dt2()
    ~~~~~~~~~~~~~~~~~~~

    Returns the minimum jerk criterion for acceleration along the path S, d²S/dt², given by:

    						d²S/dt² = S*(60t - 180t² + 120t³)

    Where S is found from S = self.find_path_length()

    """
    def __init__(self,a,b,c,d,x_initial,x_break,x_final):
    	self.a = a
    	self.b = b
    	self.c = c
    	self.d = d
    	self.x_initial = x_initial
    	self.x_break = x_break
    	self.xlim = [x_initial,x_final]
    	#self.all_values = {'A': a, 'B' : b, 'C' : c, 'D' : d, 'init' : x_initial, 'break' : x_break}
    def pp_func(self,X):
    	import numpy as np
    	result = np.piecewise(X,[X <= self.x_break, X > self.x_break], \
    		[lambda X: self.a[0] + self.b[0,0,0]*(X-self.x_initial) + self.c[0,0]*(X-self.x_initial)**2 + self.d[0,0,0]*(X-self.x_initial)**3, \
    		lambda X: self.a[1] + self.b[1,0,0]*(X-self.x_break) + self.c[1,0]*(X-self.x_break)**2 + self.d[1,0,0]*(X-self.x_break)**3])
    	return(result)
    def pp_deriv(self,X):
    	import numpy as np
    	result = np.piecewise(X,[X <= self.x_break, X > self.x_break], \
    		[lambda X: self.b[0,0,0] + 2*self.c[0,0]*(X-self.x_initial) + 3*self.d[0,0,0]*(X-self.x_initial)**2, \
    		lambda X: self.b[1,0,0] + 2*self.c[1,0]*(X-self.x_break) + 3*self.d[1,0,0]*(X-self.x_break)**2])
    	return(result)
    def pp_2deriv(self,X):
    	import numpy as np
    	result = np.piecewise(X,[X <= self.x_break, X > self.x_break], \
    		[lambda X: 2*self.c[0,0] + 6*self.d[0,0,0]*(X-self.x_initial), \
    		lambda X: 2*self.c[1,0] + 6*self.d[1,0,0]*(X-self.x_break)])
    	return(result)
    def find_max_and_min(self,x_min,x_max,y_min,y_max):
    	def find_extrema():
    		import numpy as np
    		if (4*self.c[0,0]**2 - 12*self.b[0,0,0]*self.d[0,0,0]) >= 0:
    			extrema_1 = self.x_initial + (- 2*self.c[0,0] + (4*self.c[0,0]**2 - 12*self.b[0,0,0]*self.d[0,0,0])**.5)/(6*self.d[0,0,0])
    			extrema_2 = self.x_initial + (- 2*self.c[0,0] - (4*self.c[0,0]**2 - 12*self.b[0,0,0]*self.d[0,0,0])**.5)/(6*self.d[0,0,0])
    		else:
    			extrema_1, extrema_2 = None, None
    		if (4*self.c[1,0]**2 - 12*self.b[1,0,0]*self.d[1,0,0]) >= 0:
    			extrema_3 = self.x_break + (- 2*self.c[1,0] + (4*self.c[1,0]**2 - 12*self.b[1,0,0]*self.d[1,0,0])**.5)/(6*self.d[1,0,0])
    			extrema_4 = self.x_break + (- 2*self.c[1,0] - np.sqrt(4*self.c[1,0]**2 - 12*self.b[1,0,0]*self.d[1,0,0]))/(6*self.d[1,0,0])
    		else:
    			extrema_3, extrema_4 = None, None
    		return(extrema_1,extrema_2,extrema_3,extrema_4)
    	def is_real(x_value):
    		result = not isinstance(x_value,complex)
    		return(result)
    	def is_in_appropriate_domain(x_value,x_min,x_max,segment_number):
    		if segment_number == 1:
    			result = x_value >= x_min and x_value <= self.x_break
    		elif segment_number == 2:
    			result = x_value >= self.x_break and x_value <= x_max
    		return(result)
    	def second_deriv_is_neg(x_value,segment_number):
    		if segment_number == 1:
    			x_not = self.x_initial
    		elif segment_number == 2:
    			x_not = self.x_break
    		second_deriv =  2*self.c[segment_number-1] + 6*self.d[segment_number-1]*(x_value-x_not)
    		result = second_deriv<0
    		return(result)
    	def second_deriv_is_pos(x_value,segment_number):
    		if segment_number == 1:
    			x_not = self.x_initial
    		elif segment_number == 2:
    			x_not = self.x_break
    		second_deriv =  2*self.c[segment_number-1] + 6*self.d[segment_number-1]*(x_value-x_not)
    		result = second_deriv>0
    		return(result)
    	def determine_if_max_or_min(extrema_1,extrema_2,extrema_3,extrema_4,x_min,x_max):
    		import numpy as np
    		maxima = []
    		minima = []
    		if extrema_1 != None and is_in_appropriate_domain(extrema_1,x_min,x_max,1):
    			if second_deriv_is_neg(extrema_1,1):
    				maxima.append(np.float(self.pp_func(extrema_1)))
    			elif second_deriv_is_pos(extrema_1,1):
    				minima.append(np.float(self.pp_func(extrema_1)))
    		if extrema_2 != None and is_in_appropriate_domain(extrema_2,x_min,x_max,1):
    			if second_deriv_is_neg(extrema_2,1):
    				maxima.append(np.float(self.pp_func(extrema_2)))
    			elif second_deriv_is_pos(extrema_2,1):
    				minima.append(np.float(self.pp_func(extrema_2)))
    		if extrema_3 != None and is_in_appropriate_domain(extrema_3,x_min,x_max,2):
    			if second_deriv_is_neg(extrema_3,2):
    				maxima.append(np.float(self.pp_func(extrema_3)))
    			elif second_deriv_is_pos(extrema_3,2):
    				minima.append(np.float(self.pp_func(extrema_3)))
    		if extrema_4 != None and is_in_appropriate_domain(extrema_4,x_min,x_max,2):
    			if second_deriv_is_neg(extrema_4,2):
    				maxima.append(np.float(self.pp_func(extrema_4)))
    			elif second_deriv_is_pos(extrema_4,2):
    				minima.append(np.float(self.pp_func(extrema_4)))
    		return(maxima,minima)
    	extrema_1,extrema_2,extrema_3,extrema_4 = find_extrema()
    	maxima, minima = determine_if_max_or_min(extrema_1,extrema_2,extrema_3,extrema_4,x_min,x_max)
    	return(maxima,minima)
    def is_initial_slope_positive(self,X,cutoff):
    	result = min(self.pp_deriv(X[:cutoff]))>=0
    	return(result)
    def is_within_bounds(self,x_min,x_max,y_min,y_max):
    	import numpy as np
    	maxima,minima = self.find_max_and_min(x_min,x_max,y_min,y_max)
    	if len(maxima) == 0:
    		maxima = y_max
    	if len(minima) == 0:
    		minima = y_min
    	result = np.max(maxima) <= y_max and np.min(minima) >= y_min
    	return(result)
    def print_func(self):
    	from sympy import Symbol,Lambda,pprint
    	x = Symbol('x')
    	func_1 = Lambda(x,self.a[0] + self.b[0,0,0]*(x-self.x_initial) + self.c[0,0]*(x-self.x_initial)**2 + self.d[0,0,0]*(x-self.x_initial)**3)
    	print('Function 1:\n')
    	pprint(func_1)
    	func_2 = Lambda(x,self.a[1] + self.b[1,0,0]*(x-self.x_break) + self.c[1,0]*(x-self.x_break)**2 + self.d[1,0,0]*(x-self.x_break)**3)
    	print('Function 2:\n')
    	pprint(func_2)
    def return_parameterized_X(self,t_end=1):
    	import scipy.integrate as integrate
    	import numpy as np
    	N = 1000
    	t = np.linspace(0,t_end, N + 1)
    	dt = t[1]-t[0]
    	def ode_func(x,t,t_end):
    		return(self.dS_dt(t,t_end=t_end)/np.sqrt(1 + self.pp_deriv(x)**2))
    	X = integrate.odeint(lambda x,t: ode_func(x,t,t_end),self.xlim[0],t).flatten()
    	Y = np.array(list(map(lambda x: self.pp_func(x),X)))
    	dS = np.array(list(map(lambda dx,dy: np.sqrt(dx**2+dy**2),\
    							np.gradient(X)/dt,np.gradient(Y)/dt)))
    	assert sum(abs(self.dS_dt(t,t_end=t_end)-dS))/len(dS)<1e-4, "Error in parameterizing path to dS/dt. Check ODE func."
    	return(X,Y)
    def return_parameterized_dX(self,x,t_end=1):
    	import numpy as np
    	N = 1000
    	t = np.linspace(0,t_end,N + 1)
    	# dt = t[1]-t[0]
    	dS_dt = self.dS_dt(t,t_end=t_end)
    	df_dx = self.pp_deriv(x)
    	dx_dt = np.array(list(map(lambda dS_dt,df_dx: dS_dt/np.sqrt(1 + df_dx**2),dS_dt,df_dx)))
    	dy_dt = df_dx*dx_dt
    	return(dx_dt,dy_dt)
    def return_parameterized_d2X(self,x,dx_dt,t_end=1):
    	import numpy as np
    	N = 1000
    	t = np.linspace(0,t_end,N + 1)
    	# dt = t[1]-t[0]

    	dS_dt = self.dS_dt(t,t_end=t_end)
    	d2S_dt2 = self.d2S_dt2(t,t_end=t_end)

    	df_dx = self.pp_deriv(x)
    	d2f_dx2 = self.pp_2deriv(x)

    	d2x_dt2 = np.array(list(map(lambda dx_dt,d2S_dt2,df_dx,d2f_dx2: \
    						(d2S_dt2*np.sqrt(1+df_dx**2) - df_dx*d2f_dx2*(dx_dt**2))\
    								/(1 + df_dx**2), \
    							dx_dt,d2S_dt2,df_dx,d2f_dx2)))
    	d2y_dt2 = np.array(list(map(lambda d2f_dx2,dx_dt,df_dx,d2x_dt2:  \
    						d2f_dx2*(dx_dt**2) + df_dx*(d2x_dt2),\
    							d2f_dx2,dx_dt,df_dx,d2x_dt2)))
    	return(d2x_dt2,d2y_dt2)
    def find_path_length(self):
    	import numpy as np
    	import scipy.integrate as integrate
    	return(integrate.quad(lambda x: \
    			np.sqrt(1 + self.pp_deriv(x)**2),self.xlim[0],self.xlim[1])[0])
    def dS_dt(self,t,t_end=1):
    	S_initial = 0
    	S_final = self.find_path_length()
    	return((S_final-S_initial)*(30*(t/t_end)**2 - 60*(t/t_end)**3 + 30*(t/t_end)**4)/t_end)
    def d2S_dt2(self,t,t_end=1):
    	S_initial = 0
    	S_final = self.find_path_length()
    	return((S_final-S_initial)*(60*(t/t_end) - 180*(t/t_end)**2 + 120*(t/t_end)**3)/(t_end**2))

def generate_default_path(X,Y,boundaryConditions,**kwargs):
    assert type(X)==list and len(X) == 3, "X must be a list of length 3."
    assert type(Y)==list and len(Y) == 3, "Y must be a list of length 3."
    assert type(boundaryConditions)==list and len(boundaryConditions) == 2, "boundaryConditions must be a list of length 3."
    def spline_coefficients(x1,x2,x3,y1,y2,y3,initial_slope,final_slope):
    	"""
    	Uses the values of (x1,y1), (x2,y2), and (x3,y3) to find the coefficients for the piecewise polynomial
    	equation y = a + b*(x-x_o) + c*(x-x_o)**2 + d*(x-x_o)**3) for a clamped cubic spline with one break only.
    	Returns coefficient arrays A, B, C,and D.
    	"""
    	def c_matrix(x1,x2,x3):
    		"""
    		Takes in the values of x1, x2, and x3 to create the C matrix needed to find the coefficients of a clamped
    		cubic spline with only one break (i.e. Cx = y, where x is an array of c coefficients for the
    		piecewise polynomial equation y = a + b*(x-x_o) + c*(x-x_o)**2 + d*(x-x_o)**3). Returns a matrix.
    		"""
    		import numpy as np
    		C = np.array([	[	2*(x2-x1), 		(x2-x1), 			0			],   \
    						[	(x2-x1), 		2*(x3-x1), 		(x3-x2)		],   \
    						[	0,				(x3-x2),		2*(x3-x2)	] 	], \
    						float)
    		return(C)
    	def y_vector(x1,x2,x3,y1,y2,y3,initial_slope,final_slope):
    		"""
    		Takes in the values of (x1,y1), (x2,y2), and (x3,y3) to create the y array necessary for the clamped cubic
    		spline matrix manipulation for one break only (i.e. Cx = y, where x is an array of c coefficients for the
    		piecewise polynomial equation y = a + b*(x-x_o) + c*(x-x_o)**2 + d*(x-x_o)**3). Returns an array.
    		"""
    		y = np.array([	3*(y2-y1)/(x2-x1) - 3*initial_slope ,  	\
    						3*(y3-y2)/(x3-x2) - 3*(y2-y1)/(x2-x1),  \
    						3*final_slope - 3*(y3-y2)/(x3-x2)	],  \
    						float)
    		return(y)
    	def c_coefficients(x1,x2,x3,y1,y2,y3,initial_slope,final_slope):
    		"""
    		Using matrix manipulations the equation Cx = y necessary for the c coefficients for a clamped cubic spline
    		with only one break (i.e. Cx = y, where x is an array of c coefficients for the piecewise polynomial
    		equation y = a + b*(x-x_o) + c*(x-x_o)**2 + d*(x-x_o)**3) can be rearranged such that x = C.T*y. The values
    		(x1,y1), (x2,y2), and (x3,y3) are the three points needed to the spline and initial_slope and final_slope
    		are the endpoint conditions. Returns an array.
    		"""
    		C = c_matrix(x1,x2,x3)
    		y = y_vector(x1,x2,x3,y1,y2,y3,initial_slope,final_slope)
    		CCoefficients = (np.matrix(C)**(-1))*(np.matrix(y).T)
    		return(CCoefficients)
    	def d_coefficients(x1,x2,x3,CCoefficients):
    		"""
    		Uses the c coefficients and the values of x1, x2, and x3 to find the d coefficients for the	piecewise
    		polynomial equation y = a + b*(x-x_o) + c*(x-x_o)**2 + d*(x-x_o)**3). CCoefficients must be an array with
    		three elements. Returns an array.
    		"""
    		DCoefficients = np.array([	(CCoefficients[1]-CCoefficients[0])/(3*(x2-x1)),  \
    									(CCoefficients[2]-CCoefficients[1])/(3*(x3-x2))	],  \
    									float)
    		return(DCoefficients)
    	def b_coefficients(x1,x2,x3,y1,y2,y3,CCoefficients,DCoefficients):
    		"""
    		Uses the c and d coefficients and the values of (x1,y1), (x2,y2), and (x3,y3) to find the b coefficients for
    		the	piecewise polynomial equation y = a + b*(x-x_o) + c*(x-x_o)**2 + d*(x-x_o)**3). CCoefficients must be an
    		array with two or more elements and DCoefficients must be an array with two elements. Returns an array.
    		"""
    		BCoefficients = np.array([	((y2-y1)/(x2-x1)-CCoefficients[0]*(x2-x1) - DCoefficients[0]*((x2-x1)**2)),  \
    									((y3-y2)/(x3-x2)-CCoefficients[1]*(x3-x2) - DCoefficients[1]*((x3-x2)**2)) 	]).astype(float)
    		return(BCoefficients)
    	def test_b_coefficients(x1,x2,x3,y1,y2,y3,CCoefficients,DCoefficients,expected_slope):
    		"""
    		Tests to make sure that the generated b coefficients match the expected slope. Uses the c and d coefficients
    		and the values of (x1,y1), (x2,y2), and (x3,y3) to find the b coefficients for the	piecewise polynomial
    		equation y = a + b*(x-x_o) + c*(x-x_o)**2 + d*(x-x_o)**3). CCoefficients must be an array with two or more
    		elements and DCoefficients must be an array with two elements. Returns TRUE if expected_slope equals b.
    		"""
    		B = b_coefficients(x1,x2,x3,y1,y2,y3,CCoefficients,DCoefficients)
    		result = abs(B[0]-expected_slope)< 0.001
    		return(result)
    		assert B[0]==expected_slope, "First b coefficient (%f) does not equal initial slope (%f)." (B[0],expected_slope)
    	def a_coefficients(y1,y2):
    		"""
    		Uses the y values of (x1,y1) and (x2,y2) to find the a coefficients for the	piecewise polynomial equation
    		y = a + b*(x-x_o) + c*(x-x_o)**2 + d*(x-x_o)**3). Returns an array.
    		"""
    		ACoefficients = np.array([	y1,    \
    									y2  ]).astype(float)
    		return(ACoefficients)
    	def test_endpoint_slope(b,c,d,x_n_minus_1,x_n,expected_slope):
    		"""
    		Takes in the cubic spline coefficients for the derivative of y = a + b*(x-x_n_minus_1) + c*(x-x_n_minus_1)**2 + d*(x-x_n_minus_1)**3
    		(y' = b + 2*c*(x-x_n_minus_1) + 3*d*(x-x_n_minus_1)**2)	for the last piecewise polynomial and tests to see if the expected slope at
    		the endpoint is equal to the actual	endpoint slope. The variable x_n_minus_1 is the initial value of the final piecewise polynomial
    		and x_n is the final data point. Returns TRUE if they are equal.

    		"""
    		actual_slope = b + 2*c*(x_n-x_n_minus_1) + 3*d*(x_n-x_n_minus_1)**2
    		result = abs(actual_slope-expected_slope)<0.001
    		return(result)
    	def test_for_discontinuity(a_n,b_n,c_n,d_n,x_n,x_n_plus_1,y_n_plus_1):
    		"""
    		Takes in the coefficients for a cubic spline polynomial y = a_n + b_n*(x-x_n) + c_n*(x-x_n)**2 + d_n*(x-x_n)**3
    		and tests to see if the final y value for this piecewise polynomial is equal to the initial y value of the next
    		piecewise polynomial (i.e. when x = x_n_plus_1). The variable x_n is the initial x value of the preceding
    		polynomial, and x_n_plus_1 is the transition value from one polynomial to the next. y_n_plus_1 is the initial y
    		value for the next piecewise polynomial.
    		"""
    		y_n_final = a_n + b_n*(x_n_plus_1-x_n) + c_n*(x_n_plus_1-x_n)**2 + d_n*(x_n_plus_1-x_n)**3
    		result = abs(y_n_final-y_n_plus_1)<0.001
    		return(result)

    	C = c_coefficients(x1,x2,x3,y1,y2,y3,initial_slope,final_slope)
    	D = d_coefficients(x1,x2,x3,C)
    	B = b_coefficients(x1,x2,x3,y1,y2,y3,C,D)
    	A = a_coefficients(y1,y2)

    	assert test_b_coefficients(x1,x2,x3,y1,y2,y3,C,D,initial_slope), "Initial slope does not match the expected value"
    	assert test_endpoint_slope(B[1,0,0],C[1,0],D[1,0,0],x2,x3,final_slope),"Problem with Endpoint Slope"
    	assert test_for_discontinuity(A[0],B[0,0,0],C[0,0],D[0,0,0],x1,x2,A[1]), "Jump Discontinuity at t = %f!" %x2


    	return(A,B,C[:2],D)

    initialX,randomX,finalX = X
    initialY,randomY,finalY = Y
    initialSlope,finalSlope = boundaryConditions
    A,B,C,D = spline_coefficients(initialX,randomX,finalX,initialY,randomY,finalY,initialSlope,finalSlope)
    path = Spline(A,B,C,D,initialX,randomX,finalX)

    yBounds = kwargs.get("yBounds",None)
    assert (yBounds is None) or (type(yBounds)==list and len(yBounds)==2), "yBounds must either be None or a list of length 2."
    if yBounds is not None:
        ymin = yBounds[0]
        ymax = yBounds[1]

    if yBounds is None:
        return(path)
    elif path.is_within_bounds(initialX,finalX, ymin, ymax):
        return(path)
    else:
        return(None)

def return_N_splines(N,Time):
    endpointErrorSigma = 0.005*F_MAX1
    endpointErrorTheta = 30*(np.pi/180)
    maximumDeviationInTension = 0.1*F_MAX1
    splitTimeArray = split_time_array(N,Time)

    maximumTension = maximumDeviationInTension
    minimumTension = -maximumDeviationInTension


    ValidPath = False
    count = 0
    while ValidPath == False:
        nullspaceSplineList = []
        nullspaceSplinePathList = []

        splineTimeBreaks = [0]
        splineTensionBreaks = [0]
        splineSlopeBreaks = [0]
        for i in range(N):
            splineTimeBreaks.append(splitTimeArray[i][-1])
            splineTensionBreaks.append(np.random.normal(0,endpointErrorSigma))
            if list(np.sign(splineTensionBreaks[-2:])) in [[1,1],[-1,-1]]:
                splineTensionBreaks[-1]=-splineTensionBreaks[-1]
            # splineTensionBreaks.append(0)

            splineSlopeBreaks.append(np.random.uniform(
                                        -np.tan(endpointErrorTheta),
                                        np.tan(endpointErrorTheta))
                                    )
            # splineSlopeBreaks.append(0)

        statusbar = dsb(0,N,title=return_N_splines.__name__)
        for i in range(N):
            initialTime = splineTimeBreaks[i]
            finalTime = splineTimeBreaks[i+1]

            initialTension = splineTensionBreaks[i]
            finalTension = splineTensionBreaks[i+1]

            initialSlope = splineSlopeBreaks[i]
            finalSlope = splineSlopeBreaks[i+1]

            randomTimeInBounds = False
            while randomTimeInBounds==False:
                randomTime = np.random.normal(
                        (initialTime+finalTime)/2,
                        abs(initialTime+finalTime)/8)
                if initialTime<=randomTime<=finalTime:
                    randomTimeInBounds = True
            randomTension = np.random.normal(0,endpointErrorSigma)

            path = generate_default_path(
                    [initialTime,randomTime,finalTime],\
                    [initialTension,randomTension,finalTension],\
                    [initialSlope,finalSlope],\
                    yBounds=[minimumTension,maximumTension]
            )
            if path is None:
                break
            else:
                nullspaceSplineList.append(path.pp_func(splitTimeArray[i]))
                nullspaceSplinePathList.append(path)
                del(path)
                statusbar.update(i)
        if len(nullspaceSplineList)!=0:
            nullspaceSplineArray = np.concatenate(nullspaceSplineList)
            if len(nullspaceSplineArray)==len(Time):
                ValidPath = True
            else:
                count +=1
                assert count < 100, "Trouble finding valid path after 100 trials."
    return(nullspaceSplineArray)

def return_particular_tension(t,**kwargs):
    solutionType = kwargs.get("solutionType",None)
    if solutionType is None:
        T = np.zeros(len(t))
    elif solutionType == "whiteNoise":
        T = np.random.normal(
                            loc=0.0,
                            scale=0.001*F_MAX2,
                            size=(len(t),)
                            )
    elif solutionType == "splines":
        T =return_N_splines(4,t)
    return(T)

def return_U_random_tensions(i,t,X,U,particularSolution1,**kwargs):
    """
    Takes in time scalar (float) (t), state numpy.ndarray (X) of shape (2,), and previous input numpy.ndarray (U) of shape (2,) and returns the input for this time step.

    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    **kwargs
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    1) Noise - must be an numpy.ndarray of shape (2,). Default is np.zeros((1,2)).

    2) Seed - must be a scalar value. Default is None.

    3) Bounds - must be a (2,2) list with each row in ascending order. Default is given by Tension_Bounds.

    4) MaxStep - must be a scalar (int or float). Default is MaxStep_Tension.

    """
    import random
    import numpy as np

    assert np.shape(X) == (2,) and str(type(X)) == "<class 'numpy.ndarray'>", "X must be a (2,) numpy.ndarray"
    assert np.shape(U) == (2,) and str(type(U)) == "<class 'numpy.ndarray'>", "U must be a (2,) numpy.ndarray"

    dt = t[1]-t[0]

    Noise = kwargs.get("Noise",np.zeros((2,)))
    assert np.shape(Noise) == (2,) and str(type(Noise)) == "<class 'numpy.ndarray'>", "Noise must be a (2,) numpy.ndarray"

    Seed = kwargs.get("Seed",None)
    assert type(Seed) in [float,int] or Seed is None, "Seed must be a float or an int or None."
    np.random.seed(Seed)

    Bounds = kwargs.get("Bounds",Tension_Bounds)
    assert type(Bounds) == list and np.shape(Bounds) == (2,2), "Bounds for Tension Control must be a (2,2) list."
    assert Bounds[0][0]<Bounds[0][1],"Each set of bounds must be in ascending order."
    assert Bounds[1][0]<Bounds[1][1],"Each set of bounds must be in ascending order."

    MaxStep = kwargs.get("MaxStep",MaxStep_Tension)
    assert type(MaxStep) in [int,float], "MaxStep for Tension Controller should be an int or float."

    """
    Note on QP from CVXOPT:
    Takes the form:
    min. f(x) = (1/2)*(x^T)*Q*x + p*x
    s.t. G*x<=h and A*x=b
    When calculating Q, make sure that you DO NOT INCLUDE THE 1/2! The algorithm automatically assumes this value. Therefore, always rearrange f(x) into the most reduced form shown above and then extract Q without the 1/2.
    To show functionality, run this script:
    ```py3
    import numpy as np
    import matplotlib.pyplot as plt
    import cvxopt
    cvxopt.solvers.options['show_progress'] = False
    Bounds = [[0,1],[0,1]]
    Coefficient1,Coefficient2,Constraint1 = 1,1,1
    np.random.seed()
    U = np.random.rand(2)
    I2 = cvxopt.matrix([1,0,0,1],(2,2),tc='d')
    Q1 = I2
    p1 = cvxopt.matrix([-U[0],-U[1]],tc='d')
    G1 = cvxopt.matrix([-I2,I2])
    h1 = cvxopt.matrix([-Bounds[0][0],-Bounds[1][0],Bounds[0][1],Bounds[1][1]],tc='d')
    A1 = cvxopt.matrix([Coefficient1,Coefficient2], (1,2),tc='d')
    b1 = cvxopt.matrix(Constraint1,tc='d')
    sol1 = cvxopt.solvers.qp(Q1,p1,G1,h1,A1,b1)
    X1 = [sol1['x'][0],sol1['x'][1]]
    print("Next x from min. (0.5)((x-u)^T)I(x-u):\n")
    print(X1)
    print("\n")
    Q2 = I2
    p2 = cvxopt.matrix([0,0],tc='d')
    G2 = cvxopt.matrix([-I2,I2])
    h2 = cvxopt.matrix([U[0]-Bounds[0][0],U[1]-Bounds[1][0],Bounds[0][1]-U[0],Bounds[1][1]-U[1]],tc='d')
    A2 = cvxopt.matrix([Coefficient1,Coefficient2], (1,2),tc='d')
    b2 = cvxopt.matrix(Constraint1-Coefficient1*U[0]-Coefficient2*U[1],tc='d')
    sol2 = cvxopt.solvers.qp(Q2,p2,G2,h2,A2,b2)
    X2 = [sol2['x'][0]+U[0],sol2['x'][1]+U[1]]
    print("Next x from min. (0.5)(\\xi^T)I(\\xi) where \\xi = x-u:\n")
    print(X2)
    x = np.linspace(0,1,1001)
    y = np.linspace(0,1,1001)
    cy = -Coefficient1*x/Coefficient2 + Constraint1/Coefficient2
    cx = x
    plt.plot([0,1,1,0,0],[0,0,1,1,0],'0.70',linestyle='--')
    plt.scatter([U[0]],[U[1]],c='r',marker='o')
    plt.plot(cx,cy,'b')
    plt.plot([U[0],Coefficient1*(Constraint1 - Coefficient1*U[0]-Coefficient2*U[1])/(Coefficient1**2+Coefficient2**2)+U[0]],[U[1],Coefficient2*(Constraint1 - Coefficient1*U[0]-Coefficient2*U[1])/(Coefficient1**2+Coefficient2**2)+U[1]],'r')
    plt.gca().set_aspect('equal')
    plt.scatter([Coefficient1*(Constraint1 - Coefficient1*U[0]-Coefficient2*U[1])/(Coefficient1**2+Coefficient2**2)+U[0]],[Coefficient2*(Constraint1 - Coefficient1*U[0]-Coefficient2*U[1])/(Coefficient1**2+Coefficient2**2)+U[1]],c='r',marker='o')
    plt.scatter([X1[0]],[X1[1]],c='k',marker='o')
    plt.plot([U[0],X1[0]],[U[1],X1[1]],'k')
    plt.scatter([X2[0]],[X2[1]],c='g',marker='o')
    plt.plot([U[0],X2[0]],[U[1],X2[1]],'g')
    plt.show()
    ```
    ```py3
    cvxopt.solvers.options['show_progress'] = False
    Bounds = [[0,1],[0,1]]
    Coefficient1,Coefficient2,Constraint1 = 1,1,1
    np.random.seed()
    U = np.random.rand(2)
    I2 = cvxopt.matrix([1,0,0,1],(2,2),tc='d')
    costWeights1 = [0,1]
    Q1 = sum(costWeights1)*I2
    p1 = costWeights1[1]*cvxopt.matrix([-U[0],-U[1]],tc='d')
    G1 = cvxopt.matrix([-I2,I2])
    h1 = cvxopt.matrix([-Bounds[0][0],-Bounds[1][0],Bounds[0][1],Bounds[1][1]],tc='d')
    A1 = cvxopt.matrix([Coefficient1,Coefficient2], (1,2),tc='d')
    b1 = cvxopt.matrix(Constraint1,tc='d')
    sol1 = cvxopt.solvers.qp(Q1,p1,G1,h1,A1,b1)
    X1 = [sol1['x'][0],sol1['x'][1]]
    print("Nearest Next Input:\n")
    print(X1)
    print("\n")
    costWeights2 = [1,0]
    Q2 = sum(costWeights2)*I2
    p2 = costWeights2[1]*cvxopt.matrix([-U[0],-U[1]],tc='d')
    G2 = cvxopt.matrix([-I2,I2])
    h2 = cvxopt.matrix([-Bounds[0][0],-Bounds[1][0],Bounds[0][1],Bounds[1][1]],tc='d')
    A2 = cvxopt.matrix([Coefficient1,Coefficient2], (1,2),tc='d')
    b2 = cvxopt.matrix(Constraint1,tc='d')
    sol2 = cvxopt.solvers.qp(Q2,p2,G2,h2,A2,b2)
    X2 = [sol2['x'][0],sol2['x'][1]]
    print("Minimum Total Input:\n")
    print(X2)
    print("\n")
    costWeights3 = [1,1]
    Q3 = sum(costWeights3)*I2
    p3 = costWeights3[1]*cvxopt.matrix([-U[0],-U[1]],tc='d')
    G3 = cvxopt.matrix([-I2,I2])
    h3 = cvxopt.matrix([-Bounds[0][0],-Bounds[1][0],Bounds[0][1],Bounds[1][1]],tc='d')
    A3 = cvxopt.matrix([Coefficient1,Coefficient2], (1,2),tc='d')
    b3 = cvxopt.matrix(Constraint1,tc='d')
    sol3 = cvxopt.solvers.qp(Q3,p3,G3,h3,A3,b3)
    X3 = [sol3['x'][0],sol3['x'][1]]
    print("Mixed Cost:\n")
    print(X3)
    x = np.linspace(0,1,1001)
    y = np.linspace(0,1,1001)
    cy = -Coefficient1*x/Coefficient2 + Constraint1/Coefficient2
    cx = x
    plt.plot([0,1,1,0,0],[0,0,1,1,0],'0.70',linestyle='--')
    plt.scatter([U[0]],[U[1]],c='r',marker='o')
    plt.plot(cx,cy,'b')
    plt.plot([U[0],Coefficient1*(Constraint1 - Coefficient1*U[0]-Coefficient2*U[1])/(Coefficient1**2+Coefficient2**2)+U[0]],[U[1],Coefficient2*(Constraint1 - Coefficient1*U[0]-Coefficient2*U[1])/(Coefficient1**2+Coefficient2**2)+U[1]],'r')
    plt.gca().set_aspect('equal')
    plt.scatter([Coefficient1*(Constraint1 - Coefficient1*U[0]-Coefficient2*U[1])/(Coefficient1**2+Coefficient2**2)+U[0]],[Coefficient2*(Constraint1 - Coefficient1*U[0]-Coefficient2*U[1])/(Coefficient1**2+Coefficient2**2)+U[1]],c='r',marker='o')
    plt.scatter([X1[0]],[X1[1]],c='k',marker='o')
    plt.plot([U[0],X1[0]],[U[1],X1[1]],'k')
    plt.scatter([X2[0]],[X2[1]],c='g',marker='o')
    plt.plot([U[0],X2[0]],[U[1],X2[1]],'g')
    plt.scatter([X3[0]],[X3[1]],c='b',marker='o')
    plt.plot([U[0],X3[0]],[U[1],X3[1]],'b')
    plt.show()
    ```
    Setting the weights equal to each other produces a point that is exactly inbetween the nearest next input and the point (0.5,0.5) (i.e., minimum overall activation) along the constraint line.
    """
    Coefficient1,Coefficient2,Constraint1 = return_constraint_variables(t[i],X)
    AllowableBounds_x = np.array([U[0]-MaxStep,U[0]+MaxStep])
    AllowableBounds_y = np.array([U[1]-MaxStep,U[1]+MaxStep])

    cvxopt.solvers.options['show_progress'] = False
    costWeights = [0,1] # sums to 1
    I2 = cvxopt.matrix([1.0,0.0,0.0,1.0],(2,2))
    Q = sum(costWeights)*I2
    p = costWeights[1]*cvxopt.matrix([-U[0],-U[1]],tc='d')
    # G = cvxopt.matrix([-I2,I2,-I2,I2])
    # h = cvxopt.matrix([-Bounds[0][0],-Bounds[1][0],Bounds[0][1],Bounds[1][1],\
    #             -U[0]+MaxStep,-U[1]+MaxStep,U[0]+MaxStep,U[1]+MaxStep])
    G = cvxopt.matrix([-I2,I2,-I2,I2])
    h = cvxopt.matrix([-Bounds[0][0],-Bounds[1][0],Bounds[0][1],Bounds[1][1],
                        -AllowableBounds_x[0],-AllowableBounds_y[0],
                        AllowableBounds_x[1],AllowableBounds_y[1]],tc='d')
    A = cvxopt.matrix([Coefficient1,Coefficient2], (1,2),tc='d')
    b = cvxopt.matrix(Constraint1,tc='d')
    sol=cvxopt.solvers.qp(Q, p, G, h, A, b)
    assert sol['status']=='optimal', "CVXOPT solution is not optimal."
    homogeneousSolution = np.array([sol['x'][0],sol['x'][1]])
    particularSolution2 = -R1(X)*particularSolution1/R2(X)
    particularSolution = np.array([particularSolution1,particularSolution2])
    nextU = homogeneousSolution + particularSolution
    #
    # if Constraint1 != 0:
    # 	assert Coefficient1!=0 and Coefficient2!=0, "Error with Coefficients. Shouldn't be zero with nonzero constraint."
    # else:
    # 	assert Coefficient1!=0 and Coefficient2!=0, "Error with Constraint. 0 = 0 implies all inputs valid."
    #
    #
    #
    # if Coefficient1 == 0:
    #     LowerBound_x = max(Bounds[0][0],AllowbaleBounds_x[0])
    #     UpperBound_x = min(Bounds[0][1],AllowbaleBounds_x[1])
    #     FeasibleInput1 = (UpperBound_x-LowerBound_x)*np.random.rand() + LowerBound_x
    #     FeasibleInput2 = Constraint1/Coefficient2
    # elif Coefficient2 == 0:
    #     LowerBound_y = max(Bounds[1][0],AllowableBounds_y[0])
    #     UpperBound_y = min(Bounds[1][1],AllowableBounds_y[1])
    #     FeasibleInput1 = Constraint1/Coefficient1
    #     FeasibleInput2 = (UpperBound_y-LowerBound_y)*np.random.rand() + LowerBound_y
    # else:
    #     SortedAllowableBounds = np.sort([\
    #     							(Constraint1-Coefficient2*AllowableBounds_y[0])/Coefficient1,\
    #     							(Constraint1-Coefficient2*AllowableBounds_y[1])/Coefficient1\
    #     							])
    #     SortedBounds = np.sort([(Constraint1-Coefficient2*Bounds[1][0])/Coefficient1,\
    #     							(Constraint1-Coefficient2*Bounds[1][1])/Coefficient1])
    #     LowerBound_x = max(	Bounds[0][0],\
    #      					SortedBounds[0],\
    #     					AllowableBounds_x[0],\
    #     					SortedAllowableBounds[0]\
    #     				)
    #     UpperBound_x = min(	Bounds[0][1],\
    #      					SortedBounds[1],\
    #     					AllowableBounds_x[1],\
    #     					SortedAllowableBounds[1]\
    #     				)
    #     # if UpperBound_x < LowerBound_x: import ipdb; ipdb.set_trace()
    #     assert UpperBound_x >= LowerBound_x, "Error generating bounds. Not feasible!"
    #     FeasibleInput1 = (UpperBound_x-LowerBound_x)*np.random.rand() + LowerBound_x
    #     FeasibleInput2 = Constraint1/Coefficient2 - (Coefficient1/Coefficient2)*FeasibleInput1
    #
    # return(np.array([FeasibleInput1,FeasibleInput2],ndmin=1))
    return(nextU,homogeneousSolution,particularSolution2)

def run_sim_rand_TT(N,**kwargs):
    """
    Runs one simulation for Tendon Tension control. (N = length of Time array)

    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    **kwargs
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    1) Bounds - must be a (2,2) list with each row in ascending order. Default is given by Tension_Bounds.

    2) InitialAngularAcceleration - must be a float or an int. Default is 0 (starting from rest).

    3) thresh - must be an int. Default is 25.

    """
    thresh = kwargs.get("thresh",25)
    assert type(thresh)==int, "thresh should be an int as it is the number of attempts the program should run before stopping."

    AnotherIteration = True
    AttemptNumber = 1

    while AnotherIteration == True:
        particularSolution1 = return_particular_tension(Time,solutionType=None)
        X = np.zeros((2,N))
        X_o,InitialTensions = find_initial_values_TT(**kwargs)
        X[:,0] = X_o
        U = np.zeros((2,N))
        U[:,0] = InitialTensions.T
        particularSolution = np.zeros((2,N))
        particularSolution[0,:] = particularSolution1
        particularSolution[1,0] = -R1(X[:,0])*particularSolution1[0]/R2(X[:,0])
        homogeneousSolution = np.zeros((2,N))
        homogeneousSolution[:,0] = U[:,0]

        AddNoise = False
        if AddNoise == True:
            np.random.seed(seed=None)
            NoiseArray = np.random.normal(loc=0.0,scale=0.2,size=(2,N))
        else:
            NoiseArray = np.zeros((2,N))

        try:
            cprint("Attempt #" + str(int(AttemptNumber)) + ":\n", 'green')
            statusbar = dsb(0,N-1,title=run_sim_rand_TT.__name__)
            for i in range(N-1):
                U[:,i+1],homogeneousSolution[:,i+1],particularSolution[1,i+1] = return_U_random_tensions(i,Time,X[:,i],U[:,i],
                                                    particularSolution1[i+1],
                                                    Noise=NoiseArray[:,i])
                X[:,i+1] = X[:,i] + dt*np.array([	dX1_dt(X[:,i]),\
                									dX2_dt(X[:,i],U=U[:,i+1])])
                statusbar.update(i)
            particularSolution[1,-1] = -R1(X[:,-1])*particularSolution[0,-1]/R2(X[:,-1])
            AnotherIteration = False
            return(X,U,homogeneousSolution,particularSolution)
        except:
            print('\n')
            print(" "*(get_terminal_width()\
            			- len("...Attempt #" + str(int(AttemptNumber)) + " Failed. "))\
            			+ colored("...Attempt #" + str(int(AttemptNumber)) + " Failed. \n",'red'))
            AttemptNumber += 1
            if AttemptNumber > thresh:
            	AnotherIteration=False
            	return(np.zeros((2,N)),np.zeros((2,N)),np.zeros((2,N)),np.zeros((2,N)))

def run_N_sim_rand_TT(**kwargs):
    NumberOfTrials = kwargs.get("NumberOfTrials",10)

    TotalX = np.zeros((NumberOfTrials,2,N))
    TotalU = np.zeros((NumberOfTrials,2,N))
    TotalParticularSolutions = np.zeros((NumberOfTrials,2,N))
    TotalHomogeneousSolutions = np.zeros((NumberOfTrials,2,N))
    TerminalWidth = get_terminal_width()

    print("\n")
    for j in range(NumberOfTrials):
        TrialTitle = (
            "          Trial "
            + str(j+1)
            + "/" +str(NumberOfTrials)
            + "          \n")
        print(
            " "*int(TerminalWidth/2 - len(TrialTitle)/2)
            + colored(TrialTitle,'white',attrs=["underline","bold"])
            )
        TotalX[j],TotalU[j],TotalHomogeneousSolutions[j],TotalParticularSolutions[j] = run_sim_rand_TT(N,**kwargs)

    i=0
    NumberOfSuccessfulTrials = NumberOfTrials
    while i < NumberOfSuccessfulTrials:
        if (TotalX[i]==np.zeros((2,np.shape(TotalX)[2]))).all():
            TotalX = np.delete(TotalX,i,0)
            TotalU = np.delete(TotalU,i,0)
            TotalParticularSolutions = np.delete(TotalParticularSolutions,i,0)
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
    return(TotalX,TotalU,TotalHomogeneousSolutions,TotalParticularSolutions)

def plot_N_sim_rand_TT(t,TotalX,TotalU,**kwargs):
	Return = kwargs.get("Return",False)
	assert type(Return) == bool, "Return should either be True or False"

	fig1 = plt.figure(figsize = (9,7))
	fig1_title = "Underdetermined Forced-Pendulum Example"
	plt.title(fig1_title,fontsize=16,color='gray')
	statusbar = dsb(0,np.shape(TotalX)[0],title=(plot_N_sim_rand_TT.__name__ + " (" + fig1_title +")"))
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
	statusbar.reset(title=(plot_N_sim_rand_TT.__name__ + " (" + fig2_title +")"))
	for j in range(np.shape(TotalX)[0]):
		plt.plot(t, (r(t)-TotalX[j,0,:])*180/np.pi,color='0.70')
		statusbar.update(j)
	plt.xlabel("Time (s)")
	plt.ylabel("Error (Deg)")

	statusbar.reset(
		title=(
			plot_N_sim_rand_TT.__name__
			+ " (Plotting States, Inputs, and Muscle Length Comparisons)"
			)
		)
	for j in range(np.shape(TotalX)[0]):
		if j == 0:
			fig3 = plot_states(t,TotalX[j],Return=True,InputString = "Tendon Tensions")
			fig4 = plot_inputs(t,TotalU[j],Return=True,InputString = "Tendon Tensions")
		else:
			fig3 = plot_states(t,TotalX[j],Return=True,InputString = "Tendon Tensions",\
									Figure=fig3)
			fig4 = plot_inputs(t,TotalU[j],Return=True,InputString = "Tendon Tensions", \
									Figure = fig4)
		statusbar.update(j)
	if Return == True:
		return([fig1,fig2,fig3,fig4])
	else:
		plt.show()
