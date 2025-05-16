import numpy as np
import pandas as pd
import math
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.ticker import AutoMinorLocator
import matplotlib.offsetbox as offsetbox
from scipy.interpolate import CubicHermiteSpline
from scipy.signal import savgol_filter

from model_generator import Beam

LOW_TOL = 1e-8
INTERPOLATION_POINTS = 1000

# calculate the fine discretized x-locations and their corresponding output
def calculate_values(beam, node_loc, node_x_disp, node_y_disp, node_rot):
    # approximating the transverse displacement as a cubic function, which makes sense since in FEM we are approximating the PDE w/ a cubic function
    # CubicHermiteSpline enforces slope at each node
    cs = CubicHermiteSpline(node_loc, node_y_disp, node_rot)

    # points for the smooth deformed beam curve
    x_fine = np.linspace(min(node_loc), max(node_loc), INTERPOLATION_POINTS)

    y_disp_fine = cs(x_fine)
    theta_fine = cs.derivative()(x_fine)

    # approximate the axial displacement as a linear function which also makes sense from the shape functions chosen for FEM
    x_disp_fine = np.interp(x_fine, node_loc, node_x_disp)

    # the cubic spline is continuous to the second derivative, which means the moment curve can 
    # be derived from cs. however shear is not continuous (has stepped look to it)
    # so apply a savgol filter to smooth it out
    # need to do this on a piecewise basis or else the vertical discontinuities at reactions/point loads
    # won't be vertical anymore
    moment_fine = cs.derivative().derivative()(x_fine) * beam.E * beam.I
    
    
    point_loads = beam.loads[beam.loads["type"]=="vforce"] 
    supports = beam.boundary_conditions

    discontinuity_locations = np.array([])
    for i in point_loads:
        discontinuity_locations = np.append(discontinuity_locations, i["startloc"])
    for i in supports:
        discontinuity_locations = np.append(discontinuity_locations, i["location"])
    # add start and end of beam
    discontinuity_locations = np.append(discontinuity_locations, 0.0)
    discontinuity_locations = np.append(discontinuity_locations, beam.L)

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
    
    return [x_fine, x_disp_fine, y_disp_fine, theta_fine, moment_fine, shear_fine]


# plot the deformed shape of the beam, boundary conditions, reactions
# plot table with displacement and rotation information at set intervals 
def plot_deformed_beam(beam, node_loc, x_fine, x_disp_fine, y_disp_fine, theta_fine, node_horizontal_force_rxn, node_vertical_force_rxn, 
                       node_moment_rxn, scale=1):
    
    plt.figure(figsize=(15, 9))

    # undeformed beam
    plt.plot(node_loc, np.zeros_like(node_loc), 'k--', label="Undeformed Shape")

    # deformed beam
    plt.plot(x_fine + x_disp_fine * scale, y_disp_fine * scale, 'b-', linewidth=2, label="Deformed Shape")

    plt.plot(0, 0, 'r-', linewidth=2, label="Reactions")
    plt.plot(0, 0, 'g-', linewidth=2, label="Loads")

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

    plt.text(0.02, 0.98, "Deformation factor: " + str(scale), ha='left', va='top', transform = ax.transAxes)

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
        bc_vert_force_rxn = np.round(node_vertical_force_rxn[temp_index] / 1000.0, 3)
        bc_hor_force_rxn = np.round(node_horizontal_force_rxn[temp_index] / 1000.0, 3)
        bc_moment_rxn = np.round(node_moment_rxn[temp_index] / 1000.0, 3)

        ax.add_patch(plt.arrow(bc_loc, -shape_scale*2, 0, shape_scale*2, 
                               width=0.03, length_includes_head=True, head_width=0.1, 
                               ec="none", color="red"))
        ax.text(bc_loc+0.05, -shape_scale*2, (str(bc_vert_force_rxn)+" kN"))

        # add horizontal reactions for pin and clamp
        if bc["type"]=="clamp" or bc["type"]=="pin":
            ax.add_patch(plt.arrow(bc_loc - shape_scale*2, 0.0, shape_scale*2, 0.0, 
                               width=0.03, length_includes_head=True, head_width=0.1, 
                               ec="none", color="red"))
            ax.text(bc_loc - shape_scale*3, 0.05, (str(bc_hor_force_rxn)+" kN"))

        # only add moment reactions for clamp supports
        if bc["type"]=="clamp":
            path = "arc3,rad=" + str(shape_scale*2)
            style = patches.ArrowStyle('simple', head_length=10, head_width=8, tail_width=1)
            ax.add_patch(patches.FancyArrowPatch((bc_loc+shape_scale, 0), (bc_loc-shape_scale, 0), 
                                                 arrowstyle=style, connectionstyle=path, color="red"))
            ax.text(bc_loc+0.05, shape_scale*0.75, (str(bc_moment_rxn)+" kN*m"))

    # ADD VISUAL FOR BEAM LOADING
    loads = beam.undiscretized_loads
    max_load = np.nanmax(np.abs(np.unique(np.concatenate((loads["startmag"], loads["endmag"])))))
    load_scale = shape_scale*2.0 / max_load
    for load in loads:
        startloc = load["startloc"]
        startmag = load["startmag"]
        endloc = load["endloc"]
        endmag = load["endmag"]
        if load["type"]=="force":
            #if the point load has an angle the tail needs to be rotated around the head by that angle
            headx = startloc
            heady = 0.0

            if(math.isnan(load["angle"]) == False):
                tailx = startloc - startmag*load_scale*math.sin(math.radians(load["angle"]))
                taily = -startmag*load_scale*math.cos(math.radians(load["angle"]))
            else:
                tailx = startloc
                taily = -startmag*load_scale

            ax.add_patch(plt.arrow(tailx, taily, headx - tailx, heady - taily, 
                               width=0.02, length_includes_head=True, head_width=0.1, 
                               ec="none", color="green"))
            ax.text(tailx+0.05, taily, (str(np.abs(startmag))+" kN"))
        elif load["type"]=="moment":
            path = "arc3,rad=" + str(shape_scale*2)
            style = patches.ArrowStyle('simple', head_length=10, head_width=8, tail_width=0.05)
            ax.add_patch(patches.FancyArrowPatch((startloc+shape_scale, 0), (startloc-shape_scale, 0), 
                                                 arrowstyle=style, connectionstyle=path, color="green"))
            ax.text(startloc+0.05, shape_scale*0.75, (str(startmag)+" kN*m"))
        elif load["type"]=="distributed":
            patch_coords = plt.Polygon([[startloc, 0],[endloc, 0],
                                        [endloc,-endmag*load_scale],[startloc,-startmag*load_scale]],
                                        color="green", alpha=0.5)
            ax.add_patch(patch_coords)

            if np.abs(startmag - 0.0) > LOW_TOL:
                ax.text(startloc+0.05, -startmag*load_scale, (str(np.abs(startmag))+" kN"))
            if np.abs(endmag - 0.0) > LOW_TOL:
                ax.text(endloc+0.05, -endmag*load_scale, (str(np.abs(endmag))+" kN"))

            

    # add table for node displacement and slope at specific intervals along beam
    # for table points can use the data from x_fine. get every 50th point for a table of 20 points
    factor = int(INTERPOLATION_POINTS/20.0)
    trimmed_node_loc = x_fine[::factor]
    trimmed_node_hor_disp = x_disp_fine[::factor]*1000.0
    trimmed_node_vert_disp = y_disp_fine[::factor]*1000.0
    trimmed_node_theta = theta_fine[::factor]*1000.0

    table_data = [
        ["{:.2f}".format(x) for x in trimmed_node_loc],
        ["{:.2f}".format(x_disp) for x_disp in trimmed_node_hor_disp],
        ["{:.2f}".format(y_disp) for y_disp in trimmed_node_vert_disp],
        ["{:.2f}".format(theta) for theta in trimmed_node_theta],
    ]
    table = plt.table(cellText=table_data,
                  rowLabels=["Location", "X-Disp (mm)", "Y-Disp (mm)", "Slope (rad*1e3)"], 
                  cellLoc='center', rowLoc='center',
                  loc='bottom', bbox=[0, -0.5, 1, 0.3])
    plt.subplots_adjust(bottom=0.35)

def plot_BMD_SFD(x_fine, moment_fine, shear_fine):

    fig, ax = plt.subplots(2, 1, figsize=(15, 9))
    
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

    # add table for node moment and shear at specific intervals along beam
    # for table points can use the data from x_fine. get every 50th point for a table of 20 points
    factor = int(INTERPOLATION_POINTS/20.0)
    trimmed_node_loc = x_fine[::factor]
    trimmed_node_moment = (moment_fine[::factor])/1000.0
    trimmed_node_shear = (shear_fine[::factor])/1000.0

    table_data = [
        ["{:.2f}".format(x) for x in trimmed_node_loc],
        ["{:.2f}".format(m) for m in trimmed_node_moment],
        ["{:.2f}".format(v) for v in trimmed_node_shear],
    ]
    table = plt.table(cellText=table_data,
                  rowLabels=["Location", "Moment (kNm)", "Shear (kN)"], 
                  cellLoc='center', rowLoc='center',
                  loc='bottom', bbox=[0, -0.7, 1, 0.4])
    plt.subplots_adjust(bottom=0.25, hspace=0.4)

# export the finely discretized results to an excel
def export_to_excel(x, y, theta, moment, shear):
    filepath = "./Results.xlsx"

    results = np.array([x, y, theta, moment, shear])
    results = np.transpose(results)
    
    df = pd.DataFrame(results, columns=["Node Location (m)", "Displacement (m)", "Rotation (rad)", 
                                        "Moment (kNm)", "Shear (kN)"])
    df.to_excel(filepath)
    
def run_output(beam, scale=1):

    # read output from FEA and put into numpy arrays
    results = pd.read_csv("../FEA Engine/Output/RESULTS.csv")

    node_locations = results["Node location (m)"].to_numpy()
    node_axial_displacements = results["Nodal axial displacement (m)"]
    node_transverse_displacements = results["Nodal transverse displacement (m)"].to_numpy()
    node_rotations = results["Nodal rotation (rad)"].to_numpy()
    node_horizontal_force_reactions = results["Horizontal force reactions (N)"]
    node_vertical_force_reactions = results["Vertical force reactions (N)"]
    node_moment_reactions = results["Moment reactions (Nm)"]

    x, x_disp, y_disp, theta, moment, shear = calculate_values(beam, node_locations, node_axial_displacements, 
                                                               node_transverse_displacements, node_rotations)
    
    plot_deformed_beam(beam, node_locations, x, x_disp, y_disp, theta, 
                       node_horizontal_force_reactions, node_vertical_force_reactions, node_moment_reactions, scale)
                       
    plot_BMD_SFD(x, moment, shear)

    #export_to_excel(x, y, theta, moment, shear)
    plt.savefig("PINN Convergence.pdf", format="pdf")

    plt.show()
    