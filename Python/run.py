from model_generator import Beam
from interpreter import run_output

# Write commands to generate the beam input file

''' INITIALIZE BEAM '''
# (Length, Young's Modulus, Moment of Inertia, area)
beam = Beam(10.0, 210000000000, 0.0005, 0.03)



''' APPLY BOUNDARY CONDITIONS '''
# (BC type, Location on beam)
beam.add_BC("pin", 0.0)
beam.add_BC("roller", 10.0)



''' APPLY LOADS '''
# (Load type, start location, start magnitude, end location, end magnitude)
# Point loads (force, moment) do not need an end location and magnitude, they are optional parameters
# Distributed loads require end location and magnitude
# Point force loads can have an angle specified. Angle is measured positive clockwise direction 
beam.add_load("distributed", startloc=3.0, startmag=-10000.0, endloc=7.0, endmag=-10000.0)
#beam.add_load("distributed", startloc=7.5, startmag=-3437.5, endloc=10.0, endmag=-3437.5)

#beam.add_load("force", startloc=4.0, startmag=-5000.0)
#beam.add_load("force", startloc=3.0, startmag=5000.0, angle=45.0)




''' DISCRETIZE BEAM '''
# Input approximate element length
beam.discretize(0.5)



''' CREATE INPUT AND RUN PROGRAM '''
beam.create_infile("INPUT.txt")
beam.run()


''' RUN OUTPUT FOR GRAPHING '''
# Optionally input a deformation factor for the beam deformed shape plot
run_output(beam, scale=100, plot=True)