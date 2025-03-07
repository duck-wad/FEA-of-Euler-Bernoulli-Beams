from model_generator import Beam

# Write commands to generate the beam input file

''' INITIALIZE BEAM '''
# (Length, Young's Modulus, Moment of Inertia)
beam = Beam(8.0, 200000000, 0.00012)



''' APPLY BOUNDARY CONDITIONS '''
# (BC type, Location on beam)
beam.add_BC("clamp", 0.0)
beam.add_BC("roller", 5.0)
beam.add_BC("pin", 8.0)



''' APPLY LOADS '''
# (Load type, start location, start magnitude, end location, end magnitude)
# Point loads (force, moment) do not need an end location and magnitude, they are optional parameters
# Distributed loads require end location and magnitude
beam.add_load("distributed", 0.0, 0.0, 5.0, -5000.0)
beam.add_load("moment", 5.0, 3000.0)



''' DISCRETIZE BEAM '''
# Input approximate element length
# If the mesh length is too small the C++ code spits out nonsense for some reason
# I have not figured out why
beam.discretize(0.5)

beam.create_infile("Case4.txt")