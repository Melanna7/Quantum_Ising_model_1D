"""*****************************************************************************
*
* Plot program for the outcomes of the simulation
*
*****************************************************************************"""

import os

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

#*******************************************************************************
# PARAMETERS OF THE SIMULATION
#
# SIDE_SEP = separation between the sides of different simulations.
#
#*******************************************************************************

SIDE_SEP = 1
SIDE_MIN = 4
SIDE_MAX = 9

sides = np.arange(SIDE_MIN, SIDE_MAX+1, SIDE_SEP, dtype='int')

#--- Contents ------------------------------------------------------------------

def fit_fun(x, a, b, c):
    y = c + a * np.exp(- x / b)
    return y

def load_data():
    """ Load data produced by analysis """

    data = {}
    for side in sides:
        # define data file path
        filename = f"side_{side}.dat"
        file_path = os.path.join("Data_simulations", filename)
        print("Loading " + file_path)
        # load data from each side file
        if os.path.isfile(file_path):
            data[side] = np.loadtxt(file_path, unpack='True')

    return data

#--- Fit procedure -------------------------------------------------------------

def alpha_hx(data):
    """ Fit and plot alpha(hx) behaviour"""

    alpha = []
    alp_e = []
    values_gp = {}

    # load points
    for side in sides:
        x, _, y, _, _, _, _ = data[side]
        index_nearest = min(range(len(x)), key=lambda i: abs(x[i] - 1))
        values_gp[side] = y[0:index_nearest]
    x = x[0:index_nearest]

    # fit alpha
    for idx in range(index_nearest):
        y = [values_gp[side][idx] for side in sides]
        parameters, covariance = curve_fit(fit_fun, sides, y)
        fit_b = parameters[1]
        std_deviation = np.sqrt(np.diag(covariance))
        fit_db = std_deviation[1]
        # print and store
        print(f"\nFit parameter for side {side} is: ")
        print(f"{fit_b} ± {fit_db}\n")
        alpha.append(fit_b)
        alp_e.append(fit_db)

    # plot alpha
    title = f"Behaviour of alpha(hx)"
    print("\nPlot " + title + "\n")
    # axis and style
    fig = plt.figure(title)
    plt.style.use('seaborn-whitegrid')
    plt.title(title)
    plt.ylabel(r'$ \alpha $')
    plt.xlabel('$ hx $')
    # points and function
    plt.errorbar(x, alpha, yerr=alp_e, fmt='.')
    # save and show
    plt.legend(loc='upper right')
    plt.savefig(os.path.join("Plots_and_fit", title + ".png"))
    plt.show()

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
    plt.xlabel('$ hx field $')
    # load and plot data in function of hx
    for side in sides:
        x, y, _, _, _, _, _ = data[side]
        plt.errorbar(x, y, fmt='.', label=f'side = {side}')
    # save and show
    plt.legend(loc='lower left')
    plt.savefig(os.path.join("Plots_and_fit", title + ".png"))
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
    plt.xlabel('$ hx field $')
    # load and plot data in function of hx
    for side in sides:
        x, _, y, _, _, _, _ = data[side]
        plt.errorbar(x, y, fmt='.', label=f'side = {side}')
    # save and show
    plt.legend(loc='upper left')
    plt.savefig(os.path.join("Plots_and_fit", title + ".png"))
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
    plt.xlabel('$ hx field $')
    # load and plot data in function of hx
    for side in sides:
        x, _, _, y, _, _, _ = data[side]
        plt.errorbar(x, y, fmt='.', label=f'side = {side}')
    # save and show
    plt.legend(loc='lower right')
    plt.savefig(os.path.join("Plots_and_fit", title + ".png"))
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
    plt.xlabel('$ hx field $')
    # load and plot data in function of hx
    for side in sides:
        x, _, _, _, y, _, _ = data[side]
        plt.errorbar(x, y, fmt='.', label=f'side = {side}')
    # save and show
    plt.legend(loc='upper right')
    plt.savefig(os.path.join("Plots_and_fit", title + ".png"))
    plt.show()

def plot_magnetization_x(data):
    """ Plot magnetization """

    title = "Plot magnetization X"
    print(title + "\n")
    # axis and style
    fig = plt.figure(title)
    plt.style.use('seaborn-whitegrid')
    plt.title(title)
    plt.ylabel(r'$ \langle M^x \rangle $')
    plt.xlabel('$ hx field $')
    # load and plot data in function of hx
    for side in sides:
        x, _, _, _, _, _, y = data[side]
        plt.errorbar(x, y, fmt='.', label=f'side = {side}')
    # save and show
    plt.legend(loc='lower right')
    plt.savefig(os.path.join("Plots_and_fit", title + ".png"))
    plt.show()

def plot_mag_scaling(data):
    """ Plot magnetization scaling """

    title = "Plot scaling magnetization"
    print(title + "\n")
    # axis and style
    fig = plt.figure(title)
    plt.style.use('seaborn-whitegrid')
    plt.title(title)
    plt.ylabel(r'$ \langle |M^z| \rangle * side^{1/8} $')
    plt.xlabel('$ (hx - 1) * side  $')
    # load and plot data in function of hx
    for side in sides:
        x, _, _, _, y, _, _ = data[side]
        y = y * np.power(side, (1/8))
        x = (x - 1) * side
        plt.errorbar(x, y, fmt='.', label=f'side = {side}')
    # save and show
    plt.legend(loc='upper right')
    plt.savefig(os.path.join("Plots_and_fit", title + ".png"))
    plt.show()


#-------------------------------------------------------------------------------

if __name__ == '__main__':

    data = load_data()

    plot_magnetization_z(data)
    plot_magnetization_x(data)

    plot_energy_gs(data)
    plot_energy_gap1(data)
    plot_energy_gap2(data)

    plot_mag_scaling(data)

    alpha_hx(data)
