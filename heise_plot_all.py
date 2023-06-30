"""*****************************************************************************
*
* Plot program for the outcomes of the simulation
*
*****************************************************************************"""

import os

import numpy as np
import matplotlib.pyplot as plt

#*******************************************************************************
# PARAMETERS OF THE SIMULATION
#
# SIDE_SEP = separation between the sides of different simulations.
#
#*******************************************************************************

SIDE_SEP = 1
SIDE_MIN = 4
SIDE_MAX = 15

sides = np.arange(SIDE_MIN, SIDE_MAX+1, SIDE_SEP, dtype='int')

#--- Contents ------------------------------------------------------------------

def load_data():
    """ Load data produced by analysis """

    data = {}
    for side in sides:
        # define data file path
        filename = f"side_{side}.dat"
        file_path = os.path.join("Data_Heise", filename)
        print("Loading " + file_path)
        # load data from each side file
        if os.path.isfile(file_path):
            data[side] = np.loadtxt(file_path, unpack='True')

    return data

#--- Plot procedures -----------------------------------------------------------

def plot_energy_gs(data):
    """ Plot energy of the ground state """

    title = "Plot ground state energy"
    print(title + "\n")
    # axis and style
    fig = plt.figure(title)
    plt.style.use('seaborn-whitegrid')
    plt.title(title)
    plt.ylabel(r'$ E_{GS} $')
    plt.xlabel('$ gx field $')
    # load and plot data in function of hx
    for side in sides:
        x, y, _, _, _, _, _, _, _ = data[side]
        plt.errorbar(x, y, fmt='+', label=f'side = {side}')
    # save and show
    plt.legend(loc='lower left')
    plt.savefig(os.path.join("Plots_Heisen", title + ".png"))
    plt.show()

def plot_energy_gap1(data):
    """ Plot first energy gap """

    title = "Plot first energy gap"
    print(title + "\n")
    # axis and style
    fig = plt.figure(title)
    plt.style.use('seaborn-whitegrid')
    plt.title(title)
    plt.ylabel(r'$ E_1 - E_{GS} $')
    plt.xlabel('$ gx field $')
    # load and plot data in function of hx
    for side in sides:
        x, _, y, _, _, _, _, _, _ = data[side]
        plt.errorbar(x, y, fmt='+', label=f'side = {side}')
    # save and show
    plt.legend(loc='upper left')
    plt.savefig(os.path.join("Plots_Heisen", title + ".png"))
    plt.show()

def plot_energy_gap2(data):
    """ Plot second energy gap """

    title = "Plot second energy gap"
    print(title + "\n")
    # axis and style
    fig = plt.figure(title)
    plt.style.use('seaborn-whitegrid')
    plt.title(title)
    plt.ylabel(r'$ E_2 - E_{GS} $')
    plt.xlabel('$ gx field $')
    # load and plot data in function of hx
    for side in sides:
        x, _, _, y, _, _, _, _, _ = data[side]
        plt.errorbar(x, y, fmt='+', label=f'side = {side}')
    # save and show
    plt.legend(loc='lower right')
    plt.savefig(os.path.join("Plots_Heisen", title + ".png"))
    plt.show()

def plot_magnetization_x(data):
    """ Plot magnetization X """

    title = "Plot magnetization X"
    print(title + "\n")
    # axis and style
    fig = plt.figure(title)
    plt.style.use('seaborn-whitegrid')
    plt.title(title)
    plt.ylabel(r'$ \langle M^x \rangle $')
    plt.xlabel('$ gx field $')
    # load and plot data in function of hx
    for side in sides:
        x, _, _, _, _, y, _, _ , _ = data[side]
        plt.errorbar(x, y, fmt='+', label=f'side = {side}')
    # save and show
    plt.legend(loc='lower right')
    plt.savefig(os.path.join("Plots_Heisen", title + ".png"))
    plt.show()

def plot_magnetization_y(data):
    """ Plot magnetization Y """

    title = "Plot magnetization Y"
    print(title + "\n")
    # axis and style
    fig = plt.figure(title)
    plt.style.use('seaborn-whitegrid')
    plt.title(title)
    plt.ylabel(r'$ \langle M^y \rangle $')
    plt.xlabel('$ gx field $')
    # load and plot data in function of hx
    for side in sides:
        x, _, _, _, _, _, y, _ , _ = data[side]
        plt.errorbar(x, y, fmt='+', label=f'side = {side}')
    # save and show
    plt.legend(loc='lower right')
    plt.savefig(os.path.join("Plots_Heisen", title + ".png"))
    plt.show()

def plot_magnetization_z(data):
    """ Plot magnetization """

    title = "Plot magnetization Z"
    print(title + "\n")
    # axis and style
    fig = plt.figure(title)
    plt.style.use('seaborn-whitegrid')
    plt.title(title)
    plt.ylabel(r'$ \langle |M^z| \rangle $')
    plt.xlabel('$ gx field $')
    # load and plot data in function of hx
    for side in sides:
        x, _, _, _, _, _, _, y, _ = data[side]
        plt.errorbar(x, y, fmt='+', label=f'side = {side}')
    # save and show
    plt.legend(loc='upper right')
    plt.savefig(os.path.join("Plots_Heisen", title + ".png"))
    plt.show()

def plot_order_parameter(data):
    """ Plot order parameter """

    title = "Plot order parameter"
    print(title + "\n")
    # axis and style
    fig = plt.figure(title)
    plt.style.use('seaborn-whitegrid')
    plt.title(title)
    plt.ylabel(r'$ O $')
    plt.xlabel('$ gx field $')
    # load and plot data in function of hx
    for side in sides:
        x, _, _, _, y, _, _, _, _ = data[side]
        plt.errorbar(x, y, fmt='+', label=f'side = {side}')
    # save and show
    plt.legend(loc='upper right')
    plt.savefig(os.path.join("Plots_Heisen", title + ".png"))
    plt.show()

#-------------------------------------------------------------------------------

if __name__ == '__main__':

    data = load_data()

    plot_energy_gs(data)
    plot_energy_gap1(data)
    plot_energy_gap2(data)

    plot_magnetization_z(data)
    plot_magnetization_x(data)
    plot_magnetization_y(data)

    plot_order_parameter(data)
