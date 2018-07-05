# NonlinearControl

To clone this repository, simply navigate to your desired folder and from the command line type,
```
git clone https://github.com/danhagen/NonlinearControl.git && cd NonlinearControl
```
(`pip` will be implemented in the future, but is not available right now).

## Tug of War with Springs
From a `python` interface type
```py
run Tug_of_War_with_Springs.py
```
or from the command line type
```
python Tug_of_War_with_Spring.py
```
in order to run the toy problem of controlling the position of a center mass by only applying "pulling" forces on side masses that are connected by linear springs. This problem utilizes integrator backstepping alongside a simple random number generator to select input values based on the previous input level. 

![Tug of War with Springs - Animation Screen Capture](./output_figures/Tug_of_War_with_Springs/Tug_of_War_with_Springs_Screen_Capture.PNG "Tug of War with Springs - Animation Screen Capture")


## Pendulum w/ Muscles (1 DOF, 2 DOA) 

![Pendulum with Two Muscles](./useful_figures/Muscle_forced_pendulum.PNG "Pendulum with Two Muscles")

These toy problems seek to use an integrator backstepping algorithm along with simple feasibility theory to control the nonlinear dynamics of a pendulum that is actuated by two physiologically realistic "muscles" (i.e., force producing agents that must adhere to the nonlinear viscoelasticity seen in human muscles). It is the hope of this project that underlying structure can be identified from the control dynamics that will help us better understand the challenges the nervous system faces when attempting to traverse the neuromechanical control manifold.

Control can be implemented at three levels (tendon tensions, muscle velocities, and muscle activations). Tendon-tension control and muscle-velocity control emerged from the equations of motion as reasonable "building blocks" when attempting to build the muscle-activation controller. These more upstream iterations were used primarily to ensure that the controller was producing physiologically feasible levels of state variables. To run any of these iterations you can run `IB_1DOF_2DOA_Random_Tendons_Tensions.py`,`IB_1DOF_2DOA_Random_Muscle_Velocities.py`, `IB_1DOF_2DOA_Random_Activations.py`, or `IB_1DOF_2DOA_Minimum_Activation_Velocity.py`. Note that the last iteration does not employ the random search algorithm, but instead will find the next input subject to the integration-backstepping constraint that is closest to the previous input. Each file can be altered to produce the plots and either show them or save them in the output_figures file.
