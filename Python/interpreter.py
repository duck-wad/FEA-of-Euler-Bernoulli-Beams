import numpy as np
import pandas as pd
import math
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.ticker import AutoMinorLocator
import matplotlib as mpl
from scipy.interpolate import CubicHermiteSpline

from model_generator import Beam


LOW_TOL = 1e-8

# plot the deformed shape of the beam, boundary conditions, reactions
# plot table with displacement and rotation information at set intervals 
def plot_deformed_beam(beam, node_loc, node_disp, node_slope, node_force_rxn, 
                       node_moment_rxn, scale=1.0, num_points=100):

    # approximating the displacement as a cubic function, which makes sense since in FEM we are approximating the PDE w/ a cubic function
    # CubicHermiteSpline enforces slope at each node
    cs = CubicHermiteSpline(node_loc, node_disp, node_slope)

    # points for the smooth deformed beam curve
    x_fine = np.linspace(min(node_loc), max(node_loc), num_points)
    y_fine = cs(x_fine)
    theta_fine = cs.derivative()(x_fine)

    plt.figure(figsize=(14, 8))

    # undeformed beam
    plt.plot(node_loc, np.zeros_like(node_loc), 'k--', label="Undeformed Shape")

    # deformed beam
    plt.plot(x_fine, y_fine * scale, 'b-', linewidth=2, label="Deformed Shape")

    # node displacement points
    plt.scatter(node_loc, node_disp, color='r', zorder=3, label="Nodes")

    ax = plt.gca()
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax.xaxis.set_minor_locator(AutoMinorLocator(2))

    plt.xlabel("Beam Length (m)")
    plt.ylabel("Deflection (m)")
    plt.title("Deformed Beam Shape")
    plt.legend()
    plt.axis("equal") 
    plt.grid()
    plt.grid(which="minor", linestyle="--")  

    # add boundary conditions and reactions onto the plot
    # scale shape size by the beam length. will be the diameter of circle, height of rec/triangle
    shape_scale = beam.L/25.0
    for bc in beam.boundary_conditions:
        bc_loc = bc["location"]

        # add the BC shape
        if bc["type"]=="roller":
            ax.add_patch(patches.Circle((bc_loc,-shape_scale/2.0), shape_scale/2.0, 
                                        ec="none"))
        elif bc["type"]=="clamp":
            ax.add_patch(patches.Rectangle((bc_loc-shape_scale/2.0, -shape_scale), 
                                           shape_scale, shape_scale, ec="none"))       
        elif bc["type"]=="pin":
            ax.add_patch(patches.RegularPolygon((bc_loc,-shape_scale/1.5), 3, 
                                                radius=shape_scale/1.5))

        # find the corresponding reaction for the BC
        temp_index = np.where(np.abs(bc_loc - node_loc) < LOW_TOL)[0][0]
        bc_force_rxn = node_force_rxn[temp_index]
        bc_moment_rxn = node_moment_rxn[temp_index]

        ax.add_patch(plt.arrow(bc_loc, -shape_scale*2, 0, shape_scale*2, 
                               width=0.02, length_includes_head=True, head_width=0.1, 
                               ec="none", color="red"))
        ax.text(bc_loc+0.05, -shape_scale*2-0.05, str(bc_force_rxn))

        # only add moment reactions for clamp supports
        if bc["type"]=="clamp":
            path = "arc3,rad=" + str(shape_scale*2)
            #style = patches.ArrowStyle("Fancy", head_length=0.1, head_width=0.1, tail_width=0.02)
            style = patches.ArrowStyle('simple', head_length=10, head_width=8, tail_width=0.05)
            ax.add_patch(patches.FancyArrowPatch((bc_loc+shape_scale, 0), (bc_loc-shape_scale, 0), 
                                                 arrowstyle=style, connectionstyle=path, color="red"))
            ax.text(bc_loc+0.05, shape_scale, str(bc_moment_rxn))

    # add table for node displacement and slope at specific intervals along beam
    # for table points can use the data from x_fine. get every 5th point for a table of 20 points
    trimmed_node_loc = x_fine[::5]
    trimmed_node_disp = y_fine[::5]
    trimmed_node_theta = theta_fine[::5]

    table_data = [
        ["{:.2f}".format(x) for x in trimmed_node_loc],
        ["{:.2f}".format(y) for y in trimmed_node_disp],
        ["{:.2f}".format(theta) for theta in trimmed_node_theta],
    ]
    table = plt.table(cellText=table_data,
                  rowLabels=["Location", "Displacement (m)", "Slope (rad)"], 
                  cellLoc='center', rowLoc='center',
                  loc='bottom', bbox=[0, -0.5, 1, 0.3])
    plt.subplots_adjust(bottom=0.35)

# take a list of location input and return the displacement, rotation etc.
# use same CubicHermitSpline to interpolate between points
def query_point(loc, node_locations, node_displacements, node_rotations):
    pass

def run_output(beam):

    # read output from FEA and put into numpy arrays
    results = pd.read_csv("../FEA Engine/Output/RESULTS.csv")

    node_locations = results["Node location (m)"].to_numpy()
    node_displacements = results["Nodal displacement (m)"].to_numpy()
    node_rotations = results["Nodal rotation (rad)"].to_numpy()
    node_force_reactions = results["Force reactions (N)"]
    node_moment_reactions = results["Moment reactions (Nm)"]

    plot_deformed_beam(beam, node_locations, node_displacements, node_rotations, 
                       node_force_reactions, node_moment_reactions, scale=1)
    
    plt.show()
