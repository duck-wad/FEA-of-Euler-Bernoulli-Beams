import numpy as np

from model_generator import Beam
from interpreter import run_output

# create the range of loads
loads = np.arange(1000, 20001, 1000)

for value in loads:

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
    beam.add_load("distributed", startloc=3.0, startmag=-1.0*value, endloc=7.0, endmag=-1.0*value)



    ''' DISCRETIZE BEAM '''
    # Input approximate element length
    beam.discretize(0.5)



    ''' CREATE INPUT AND RUN PROGRAM '''
    beam.create_infile("INPUT.txt")
    beam.run()



    ''' RUN OUTPUT FOR GRAPHING '''
    # Optionally input a deformation factor for the beam deformed shape plot
    run_output(beam, scale=100, filename=str(value), plot=False)