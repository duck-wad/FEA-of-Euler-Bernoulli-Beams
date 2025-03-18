import numpy as np
import pandas as pd
import math
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.ticker import AutoMinorLocator
import matplotlib as mpl
from scipy.interpolate import CubicHermiteSpline, CubicSpline
from scipy.signal import savgol_filter

from model_generator import Beam

LOW_TOL = 1e-8
INTERPOLATION_POINTS = 1000

# calculate the fine discretized x-locations and their corresponding output
def calculate_values(beam, node_loc, node_disp, node_rot):
    # approximating the displacement as a cubic function, which makes sense since in FEM we are approximating the PDE w/ a cubic function
    # CubicHermiteSpline enforces slope at each node
    cs = CubicHermiteSpline(node_loc, node_disp, node_rot)

    # points for the smooth deformed beam curve
    x_fine = np.linspace(min(node_loc), max(node_loc), INTERPOLATION_POINTS)

    y_fine = cs(x_fine)
    theta_fine = cs.derivative()(x_fine)

    # the cubic spline is continuous to the second derivative, which means the moment curve can 
    # be derived from cs. however shear is not continuous (has stepped look to it)
    # so apply a savgol filter to smooth it out
    # need to do this on a piecewise basis or else the vertical discontinuities at reactions/point loads
    # won't be vertical anymore
    moment_fine = cs.derivative().derivative()(x_fine) * beam.E * beam.I

    point_loads = beam.loads[beam.loads["type"]=="force"] 
    supports = beam.boundary_conditions

    discontinuity_locations = np.array([])
    for i in point_loads:
        discontinuity_locations = np.append(discontinuity_locations, i["startloc"])
    for i in supports:
        discontinuity_locations = np.append(discontinuity_locations, i["location"])

    discontinuity_locations = np.sort(np.unique(discontinuity_locations))

    # get the indices in x_fine where the discontinuities should occur
    discontinuity_indices = np.array([])
    for i in discontinuity_locations:
        exact_loc = np.where(x_fine == i)[0]
        if exact_loc.size != 0:
            discontinuity_indices = np.append(discontinuity_indices, exact_loc[0])
            continue
        filtered_x_fine = np.where(x_fine < i)[0]
        discontinuity_indices = np.append(discontinuity_indices, np.argmax(filtered_x_fine))

    # splice x_fine into chunks defined by the discontinuities
    shear_fine = np.array([])
    for i in range(len(discontinuity_indices)-1):
        # add 1 to the index so the starting location always to the right of the discontinuity node
        # except at the start of the beam, in which case the index is included in the start
        start = int(discontinuity_indices[i]+1)
        end = int(discontinuity_indices[i+1]+1)
        if i == 0:
            start = int(discontinuity_indices[i])
        x_spliced = x_fine[start:end]
        
        # apply the savgol filter on a piecewise basis to preserve vertical discontinuities
        # apply filter twice for better smoothing
        shear_spliced = savgol_filter(savgol_filter(
            cs.derivative().derivative().derivative()(x_spliced), 
            window_length=99, polyorder=3), window_length=99, polyorder=1)*beam.E*beam.I
        
        shear_fine = np.append(shear_fine, shear_spliced)   

    return [x_fine, y_fine, theta_fine, moment_fine, shear_fine]


# plot the deformed shape of the beam, boundary conditions, reactions
# plot table with displacement and rotation information at set intervals 
def plot_deformed_beam(beam, node_loc, x_fine, y_fine, theta_fine, node_force_rxn, 
                       node_moment_rxn, scale=1):
    
    plt.figure(figsize=(14, 8))

    # undeformed beam
    plt.plot(node_loc, np.zeros_like(node_loc), 'k--', label="Undeformed Shape")

    # deformed beam
    plt.plot(x_fine, y_fine * scale, 'b-', linewidth=2, label="Deformed Shape")

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

    # ADD VISUAL FOR BEAM LOADING

    # add table for node displacement and slope at specific intervals along beam
    # for table points can use the data from x_fine. get every 50th point for a table of 20 points
    factor = int(INTERPOLATION_POINTS/20.0)
    trimmed_node_loc = x_fine[::factor]
    trimmed_node_disp = y_fine[::factor]
    trimmed_node_theta = theta_fine[::factor]

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

def plot_BMD_SFD(beam, x_fine, moment_fine, shear_fine):

    fig, ax = plt.subplots(2, 1, figsize=(14, 8))
    
    # plot in kNm
    ax[0].plot(x_fine, moment_fine/1000, 'r-', linewidth=2, label="Moment")
    ax[0].plot(x_fine, np.zeros_like(x_fine), 'k--', label="Undeformed Shape")
    ax[0].grid()
    ax[0].set_title("Bending Moment Diagram")
    ax[0].set_xlabel("Beam Length (m)")
    ax[0].set_ylabel("Moment (kN*m)")
    ax[0].grid(which="minor", linestyle="--")
    ax[0].yaxis.set_minor_locator(AutoMinorLocator(2))
    ax[0].xaxis.set_minor_locator(AutoMinorLocator(2))

    # plot in kN
    ax[1].plot(x_fine, shear_fine/1000, 'g-', linewidth=2, label="Shear")
    ax[1].plot(x_fine, np.zeros_like(x_fine), 'k--', label="Undeformed Shape")
    ax[1].grid()
    ax[1].set_title("Shear Force Diagram")
    ax[1].set_xlabel("Beam Length (m)")
    ax[1].set_ylabel("Shear (kN)")
    ax[1].grid(which="minor", linestyle="--")
    ax[1].yaxis.set_minor_locator(AutoMinorLocator(2))
    ax[1].xaxis.set_minor_locator(AutoMinorLocator(2))

    # add table similar to the displacement table in deformed shape function

    fig.tight_layout()


# take a list of location input and return the displacement, rotation etc.
# use same CubicHermitSpline to interpolate between points
def query_point(loc, node_loc, node_disp, node_rot):
    
    cs = CubicHermiteSpline(node_loc, node_disp, node_rot)
    x_fine = np.linspace(min(node_loc), max(node_loc), 100)
    y_fine = cs(x_fine)
    theta_fine = cs.derivative()(x_fine)

def run_output(beam):

    # read output from FEA and put into numpy arrays
    results = pd.read_csv("../FEA Engine/Output/RESULTS.csv")

    node_locations = results["Node location (m)"].to_numpy()
    node_displacements = results["Nodal displacement (m)"].to_numpy()
    node_rotations = results["Nodal rotation (rad)"].to_numpy()
    node_force_reactions = results["Force reactions (N)"]
    node_moment_reactions = results["Moment reactions (Nm)"]

    x, y, theta, moment, shear = calculate_values(beam, node_locations, node_displacements, node_rotations)
    
    plot_deformed_beam(beam, node_locations, x, y, theta, 
                       node_force_reactions, node_moment_reactions, 10)
                       
    plot_BMD_SFD(beam, x, moment, shear)
    

    plt.show()
