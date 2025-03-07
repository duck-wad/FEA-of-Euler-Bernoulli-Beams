# Beam class
import numpy as np
import subprocess
import os

LOW_TOL = 1e-8
BCtype = [("type", "U10"), ("location", "f8")]
loadtype = [("type", "U12"), ("startloc", "f8"), ("endloc", "f8"), ("startmag", "f8"), ("endmag", "f8")]

class Beam:
    def __init__(self, L, E, I):
        self.L = L
        self.E = E
        self.I = I

        self.node_list = np.array([0.0, L])
        self.boundary_conditions = np.empty(0, dtype=BCtype)
        self.loads = np.empty(0, dtype=loadtype)

    def print_info(self):
        print("Length: " + str(self.L))
        print("Young's Modulus: " + str(self.E))
        print("Moment of Inertia: " + str(self.I))
        print("Nodes: " + str(self.node_list))
        print("Boundary Conditions: " + str(self.boundary_conditions))
        print("Applied Loads : " + str(self.loads))

    def add_BC(self, type, location):

        if type not in ("clamp", "roller", "pin"):
            print("Boundary \"" + type + "\" is not a valid boundary type.")
            return
        if (location > self.L or location < 0):
            print("Boundary is not located on beam.")
            return
        
        # check if BC has already been applied at this location
        temp_index = np.where(np.abs(self.boundary_conditions["location"] - location) < LOW_TOL)[0]
        if temp_index.size > 0:
            print("Warning: BC at the location " + str(location) + " has been detected. It will be overridden.")
            self.boundary_conditions = np.delete(self.boundary_conditions, temp_index)

        self.boundary_conditions = np.append(self.boundary_conditions, np.array([(type, location)], dtype=self.boundary_conditions.dtype))

    def add_load(self, type, startloc, startmag, endloc=None, endmag=None):

        # safety checks
        if type not in ("moment", "force", "distributed"):
            print("Load \"" + type + "\" is not a valid load type.")
            return
        if ((startloc > self.L or startloc < 0) or (endloc != None and (endloc > self.L or endloc < 0))):
            print("Load is not located on beam.")
            return
        if(startloc == endloc):
            print("Start location and end location cannot be the same.")
            return
        if (type == "moment" and endloc != None):
            print("Moment must be applied at a single point.")
            return
        if (endloc != None and endmag == None):
            print("Loads with an end location must have a specified end magnitude.")
            return
        if (type == "distributed" and endloc==None):
            print("Distributed loads must have an end location.")
            return

        self.loads = np.append(self.loads, np.array([(type, startloc, endloc, startmag, endmag)], dtype=self.loads.dtype))

    # in order to discretize, first get a list of all the points where BC or load has been applied
    # ensure those points are nodes in the node list
    # sort the node list from left to right (assuming x-positive to the right)
    # for each adjacent point pair in the node list, discretize the segment between with the approximate mesh length
    def discretize(self, meshlength):
        bc_points = self.boundary_conditions["location"]
        load_points = np.unique(np.append(self.loads["startloc"], self.loads["endloc"]))
        # clean nan values
        load_points = load_points[~np.isnan(load_points)]
        self.node_list = np.unique(np.concatenate((self.node_list, bc_points, load_points)))
        self.node_list = np.sort(self.node_list)

        discretized_list = np.zeros(0)
        for i in range(len(self.node_list)-1):
            start = self.node_list[i]
            end = self.node_list[i+1]
            diff = end - start
            numel = np.abs(round(diff / meshlength))
            segment_meshlength = diff / numel

            segment_nodes = np.zeros(0)
            increment = start
            while increment < end:
                segment_nodes = np.append(segment_nodes, increment)
                increment += segment_meshlength
            discretized_list = np.append(discretized_list, segment_nodes)
                
        self.node_list = np.unique(np.append(self.node_list, discretized_list))

        self.discretize_distributed_loads()         

    # after discretizing the beam, any distributed loads applied to the beam must also be discretized
    def discretize_distributed_loads(self):
        distributed_loads = self.loads[self.loads["type"] == "distributed"]
        # remove distributed loads from the original loads list
        self.loads = self.loads[self.loads["type"] != "distributed"]

        discretized_loads = np.empty(0, dtype=loadtype)
        # loop over the list of distributed loads
        # get the number of elements in the range of distributed load
        # divide the change in load by the number of elements to get load change per element
        for load in distributed_loads:
            
            start_node = load["startloc"]
            end_node = load["endloc"]
            start_mag = load["startmag"]
            end_mag = load["endmag"]
            if (load["startloc"] > load["endloc"]):
                start_node = load["endloc"]
                end_node = load["startloc"]
                start_mag = load["endmag"]
                end_mag = load["startmag"]

            segment_nodes = self.node_list[(self.node_list >= start_node) & (self.node_list <= end_node)]

            # number of elements in the range is 1 less than the number of nodes
            numel = len(segment_nodes)-1

            load_change = (end_mag - start_mag) / numel

            temp_startmag = start_mag
            segment_loads = np.zeros(0, dtype=loadtype)
            for i in range(numel):
                segment_loads = np.append(segment_loads, np.array([("distributed", segment_nodes[i], segment_nodes[i+1], temp_startmag, temp_startmag+load_change)], dtype=loadtype))
                temp_startmag += load_change

            discretized_loads = np.append(discretized_loads, segment_loads)
        
        self.loads = np.append(self.loads, discretized_loads)

    def create_infile(self):
        f = open("../FEA Engine/Input/INPUT.txt", "w")
        
        # write materials list
        f.write("[MATERIALS]" + "\n")
        f.write("E: " + str(self.E) + " I: " + str(self.I) + "\n \n")

        # write nodes
        f.write("[NODES]" + "\n")
        f.write("numnodes: " + str(len(self.node_list)) + "\n")
        for i in range(len(self.node_list)):
            f.write("node" + str(i+1) + " coord: " + str(self.node_list[i]) + "\n")
        f.write("\n")

        # write elements
        # elements numbered left to right along beam
        f.write("[ELEMENTS]" + "\n")
        f.write("numelems: " + str(len(self.node_list)-1) + "\n")
        for i in range(len(self.node_list)-1):
            f.write("elem" + str(i+1) + " start: " + str(i+1) + " end: " + str(i+2) + "\n")
        f.write("\n")
        
        # write boundary conditions
        f.write("[BOUNDARY_CONDITIONS]" + "\n")
        f.write("numbcs: " + str(len(self.boundary_conditions)) + "\n")
        for i in range(len(self.boundary_conditions)):
            # find the corresponding node to the location
            # since node list is sorted left to right, node is the location in the node_list array
            node = np.where(np.abs(self.node_list - self.boundary_conditions[i]["location"]) < LOW_TOL)[0][0] + 1
            
            f.write("bc" + str(i+1) + " type: " + str(self.boundary_conditions[i]["type"]) + " location: " + str(node) + "\n")
        f.write("\n")

        # write loads
        f.write("[LOADS]" + "\n")
        f.write("numloads: " + str(len(self.loads)) + "\n")
        for i in range(len(self.loads)):

            startnode = np.where(np.abs(self.node_list - self.loads[i]["startloc"]) < LOW_TOL)[0][0] + 1

            if(self.loads[i]["type"] == "distributed"):
                endnode = np.where(np.abs(self.node_list - self.loads[i]["endloc"]) < LOW_TOL)[0][0] + 1
                f.write("load" + str(i+1) + " type: " + str(self.loads[i]["type"]) + " start: " + str(startnode) + " end: " + str(endnode) + " startmagnitude: " + str(self.loads[i]["startmag"]) + " endmagnitude: " + str(self.loads[i]["endmag"]) + "\n")
            else:
                f.write("load" + str(i+1) + " type: " + str(self.loads[i]["type"]) + " location: " + str(startnode) + " magnitude: " + str(self.loads[i]["startmag"]) + "\n")
        f.close()

    def run(self):
        # user shouldn't have to rebuild the executable for every run
        '''
        build_result = subprocess.run(["bash", "../FEA Engine/build.sh"], cwd="../FEA Engine", capture_output=True, text=True, shell=True)

        print("Build Output:\n", build_result.stdout)
        if build_result.returncode != 0:
            print("Build Failed:\n", build_result.stderr)
            exit(1)
        '''

        cpp_executable = "../FEA Engine/program.exe"  

        print("Running C++ program...")
        try:
            subprocess.run([cpp_executable], cwd="../FEA Engine", check=True)
            print(f"C++ program completed.")
        except FileNotFoundError:
            print(f"Error: C++ executable '{cpp_executable}' not found.")
        except subprocess.CalledProcessError as e:
            print(f"Error: C++ program failed with error code {e.returncode}.")
