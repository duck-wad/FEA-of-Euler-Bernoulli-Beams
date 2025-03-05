# Beam class
import numpy as np

class Beam:
    def __init__(self, L, E, I):
        self.L = L
        self.E = E
        self.I = I

        self.node_list = np.array([0.0, L])
        self.boundary_conditions = np.empty(0, dtype=[("location", "f8"), ("type", "U10")])

    def print_info(self):
        print("Length: " + str(self.L))
        print("Young's Modulus: " + str(self.E))
        print("Moment of Inertia: " + str(self.I))
        print("Nodes: " + str(self.node_list))
        print("Boundary Conditions: " + str(self.boundary_conditions))

    def add_BC(self, location, type):

        if type not in ("clamp", "roller", "pin"):
            print("Boundary \"" + type + "\" is not a valid boundary type.")
            return
        if (location > self.L or location < 0):
            print("Boundary is not located on beam.")
            return
        
        # check if BC has already been applied at this location
        temp_index = np.where(self.boundary_conditions["location"] == location)[0]
        if temp_index.size > 0:
            print("Warning: BC at the location " + str(location) + " has been detected. It will be overridden.")
            self.boundary_conditions = np.delete(self.boundary_conditions, temp_index)

        self.boundary_conditions = np.append(self.boundary_conditions, np.array([(location, type)], dtype=self.boundary_conditions.dtype))


beam = Beam(8.0, 200000000, 0.00012)
beam.add_BC(0.0, "clamp")
beam.add_BC(5.0, "roller")
beam.add_BC(8.0, "pin")

beam.print_info()
    






# Generate beam boundary conditions

# Generate beam loads

# Discretize beam

