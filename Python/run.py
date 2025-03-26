from model_generator import Beam
from interpreter import run_output

# Write commands to generate the beam input file

''' INITIALIZE BEAM '''
# (Length, Young's Modulus, Moment of Inertia)
beam = Beam(8.0, 2e11, 0.00012, 0.00003)



''' APPLY BOUNDARY CONDITIONS '''
# (BC type, Location on beam)
beam.add_BC("clamp", 0.0)
beam.add_BC("roller", 5.0)
beam.add_BC("pin", 8.0)



''' APPLY LOADS '''
# (Load type, start location, start magnitude, end location, end magnitude)
# Point loads (force, moment) do not need an end location and magnitude, they are optional parameters
# Distributed loads require end location and magnitude
# Point force loads can have an angle specified. Angle is measured positive clockwise direction 
beam.add_load("distributed", startloc=0.0, startmag=0.0, endloc=5.0, endmag=-5000.0)
beam.add_load("moment", startloc=5.0, startmag=3000.0)
beam.add_load("force", startloc=4.0, startmag=-5000.0)
beam.add_load("force", startloc=7.0, startmag=5000.0, angle=45.0)



''' DISCRETIZE BEAM '''
# Input approximate element length
beam.discretize(0.1)



''' CREATE INPUT AND RUN PROGRAM '''
beam.create_infile("INPUT.txt")
beam.run()


''' RUN OUTPUT FOR GRAPHING '''
# Optionally input a deformation factor for the beam deformed shape plot
run_output(beam, 1000)