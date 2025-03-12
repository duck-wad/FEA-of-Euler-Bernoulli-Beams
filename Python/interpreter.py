import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import CubicHermiteSpline
from model_generator import Beam
from matplotlib.ticker import AutoMinorLocator

def plot_deformed_beam(node_loc, node_disp, node_slope, scale=1.0, num_points=100):

    # approximating the displacement as a cubic function, which makes sense since in FEM we are approximating the PDE w/ a cubic function
    # CubicHermiteSpline enforces slope at each node
    cs = CubicHermiteSpline(node_loc, node_disp, node_slope)

    # points for the smooth deformed beam curve
    x_fine = np.linspace(min(node_loc), max(node_loc), num_points)
    y_fine = cs(x_fine)

    plt.figure(figsize=(14, 8))

    # undeformed beam
    plt.plot(node_loc, np.zeros_like(node_loc), 'k--', label="Undeformed Shape")

    # deformed beam
    plt.plot(x_fine, y_fine * scale, 'b-', linewidth=2, label="Deformed Shape")

    # node displacement points
    #plt.scatter(node_loc, node_disp, color='r', zorder=3, label="Nodes")

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

    # add boundary condition shapes onto the plot

    # add reaction values and arrows onto the plot

    # add table for node displacement and slope at specific intervals along beam

""" def run_output(beam):
    plot_deformed_beam(node_locations, node_displacements, node_rotations, scale=1)
    plt.show()
 """

def output():

    # read output from FEA and put into numpy arrays
    results = pd.read_csv("../FEA Engine/Output/RESULTS.csv")

    node_locations = results["Node location (m)"].to_numpy()
    node_displacements = results["Nodal displacement (m)"].to_numpy()
    node_rotations = results["Nodal rotation (rad)"].to_numpy()

    plot_deformed_beam(node_locations, node_displacements, node_rotations, scale=1)
    plt.show()


output()