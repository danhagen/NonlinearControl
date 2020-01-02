import numpy as np

# params dictionary

params = {
    "Joint Inertia" : 1.15e-2, # kg⋅m²
    "Joint Damping" : 0.001, # N⋅s⋅m⁻¹
    "Joint Mass" : 0.541, # kg
    "Joint Moment Arm" : 0.05, # m
    "Link Center of Mass" : 0.085, # m
    "Link Length" : 0.3, # m
    "Motor Inertia" : 6.6e-5, # kg⋅m²
    "Motor Damping" : 0.00462, # N⋅s⋅m⁻¹
    "Motor Moment Arm" : 0.02, # m
    "Spring Stiffness Coefficient" : 200, # N
    "Spring Shape Coefficient" : 10, # unit-less
    "Simulation Duration" : 100, # s
    "dt" : 0.01 # s
}

# h is the step used to determine the derivative
h = 0.000001

#gravity
gr = 9.81 # m/s²
# gr = 0
