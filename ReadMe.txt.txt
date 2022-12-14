The submitted files include:
-The 8 figures
-A module for the Ball class and a Container class derived from the ball class
-A module for a Simulation class which is a compostion of the Ball class and the Container class
-A module for the testing of the class and for the investigation of physics called "Main.py"

To work with the simulation only the "Main.py" module needs to be used.
The first block of code which imports the functions needs to be run
The second block gives some parameters for graphs.

To create a simulation, the third block can be run, or a variable can be created as following:
Variable* = Simulation(mass*, radius*, numberofballs*, containerradius*, containermass*, sigma*)
Everything with a star should be replaced to the parameters you want.
mass gives the mass of the balls in kg, set to 1.0 kg if not given an argument
radius gives the radius of the balls in m, set to 1.0 m if not given an argument
numberofballs gives the number of the balls used in the simulation, set to 10 if not given an argument
containerradius gives the radius of the container in m, set to 10.0 m if not given an argument
containermass gives the mass of the balls in kg, set to 10^20 kg if not given an argument
sigma gives the Standard deviation of the velocity of the balls, set to 20.0 ms^-1 if not given an argument

If the error 'list is out of index' comes up, then the balls don't fit in the container.

To run a simulation a certain number of frames
Variable*.run(frames, Animate)
if Animate is true, then a animation of the collision will run 
to see this properly run the command %matplotlib auto in the console before running the command

Each simulation has many functions to investigate the various parameters.
The main file has the most relevant parameter checks, but if you want you can read into the specfic avaiable functions 
in the simulation module, the Thermodynamic ball module 
or by simply typing Class*? into the console where Class* is the class you want to find
