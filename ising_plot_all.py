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
SIDE_MAX = 14

plt.style.use('ggplot')

sides = np.arange(SIDE_MIN, SIDE_MAX+1, SIDE_SEP, dtype='int')

#--- Contents ------------------------------------------------------------------

def fit_pow(x, a, b):
    y = b * np.power(x, - a)
    return y

def fit_exp(x, a, b):
    y = b * np.exp(- a * x)
    return y

def load_data():
    """ Load data produced by analysis """

    data = {}
    for side in sides:
        # define data file path
        filename = f"side_{side}.dat"
        file_path = os.path.join("Data_Ising", filename)
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
    plt.xlabel('transverse field')
    # load and plot data in function of hx
    for side in sides:
        x, y, _, _, _, _, _, _, _ = data[side]
        plt.errorbar(x, y, fmt='.', label=f'side = {side}')
    # save and show
    plt.legend(loc='lower left')
    plt.savefig(os.path.join("Plots_and_fit", title + ".png"))
    plt.show()

def plot_energy_gap1(data):
    """ Plot first energy gap """

    title = "Plot energy first gap"
    print(title + "\n")
    # axis and style
    fig = plt.figure(title)
    plt.style.use('seaborn-whitegrid')
    plt.title(title)
    plt.ylabel(r'$ E_1 - E_{GS} $')
    plt.xlabel('transverse field')
    # load and plot data in function of hx
    for side in sides:
        x, _, y, _, _, _, _, _, _ = data[side]
        plt.errorbar(x, y, fmt='.', label=f'side = {side}')
    # save and show
    plt.legend(loc='upper left')
    plt.savefig(os.path.join("Plots_and_fit", title + ".png"))
    plt.show()

def plot_energy_gap2(data):
    """ Plot second energy gap """

    title = "Plot energy second gap"
    print(title + "\n")
    # axis and style
    fig = plt.figure(title)
    plt.style.use('seaborn-whitegrid')
    plt.title(title)
    plt.ylabel(r'$ E_2 - E_{GS} $')
    plt.xlabel('transverse field')
    # load and plot data in function of hx
    for side in sides:
        x, _, _, y, _, _, _, _, _ = data[side]
        plt.errorbar(x, y, fmt='.', label=f'side = {side}')
    # save and show
    plt.legend(loc='lower right')
    plt.savefig(os.path.join("Plots_and_fit", title + ".png"))
    plt.show()

#--- Fit procedure -------------------------------------------------------------

def alpha_hx(data):
    """ Fit and plot alpha(hx) behaviour"""

    alpha = []
    alp_e = []
    values_gp = {}

    # load points (they are stored in reverse order from 2.0 to 0.)
    x, _, _, _, _, _, _, _, _ = data[10]
    index_nearest = min(range(len(x)), key=lambda i: abs(x[i] - 1))
    print(f"Nearest index to 1 is {index_nearest} with hx = {x[index_nearest]}")

    x = x[index_nearest:]
    for side in sides:
        _, _, y, _, _, _, _, _, _ = data[side]
        values_gp[side] = y[index_nearest:]

    # fit alpha
    for idx in range(len(x)):
        y = [values_gp[side][idx] for side in sides]
        if(idx == 0):
            fit_fun = fit_pow
        else:
            fit_fun = fit_exp
        parameters, covariance = curve_fit(fit_fun, sides, y)
        fit_a = parameters[0]
        std_deviation = np.sqrt(np.diag(covariance))
        fit_da = std_deviation[0]
        # print and store
        print(f"\nFit parameter for {x[idx]} is: ")
        print(f"{fit_a} ± {fit_da}\n")
        if (idx > 0) :
            alpha.append(fit_a)
            alp_e.append(fit_da)
    x = x[1:]

    # plot alpha
    title = f"Behaviour of alpha(hx)"
    print("\nPlot " + title + "\n")
    # axis and style
    fig = plt.figure(title)
    plt.style.use('seaborn-whitegrid')
    plt.title(title)
    plt.ylabel(r'$ \alpha $')
    plt.xlabel('transverse field')
    # points and function
    plt.errorbar(x, alpha, yerr=alp_e, fmt='.')
    # save and show
    plt.savefig(os.path.join("Plots_and_fit", title + ".png"))
    plt.show()

#--- Magn procedures -----------------------------------------------------------

def plot_magnetization_z(data):
    """ Plot magnetization """

    title = "Plot magnetization Z"
    print(title + "\n")
    # axis and style
    fig = plt.figure(title)
    plt.style.use('seaborn-whitegrid')
    plt.title(title)
    plt.ylabel(r'order parameter $ M^z $')
    plt.xlabel('transverse field')
    # load and plot data in function of hx
    for side in sides:
        x, _, _, _, y, _, _, _, _ = data[side]
        plt.errorbar(x, y, fmt='.', label=f'side = {side}')
    # save and show
    plt.legend(loc='upper right')
    plt.savefig(os.path.join("Plots_and_fit", title + ".png"))
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
    plt.xlabel('transverse field')
    # load and plot data in function of hx
    for side in sides:
        x, _, _, _, _, y, _, _ , _ = data[side]
        plt.errorbar(x, y, fmt='.', label=f'side = {side}')
    # save and show
    plt.legend(loc='lower right')
    plt.savefig(os.path.join("Plots_and_fit", title + ".png"))
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
    plt.xlabel('transverse field')
    # load and plot data in function of hx
    for side in sides:
        x, _, _, _, _, _, y, _ , _ = data[side]
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
    plt.ylabel(r'$ \langle |M^z| \rangle * N^{\beta / \nu} $')
    plt.xlabel('$ (hx - 1) * side  $')
    # load and plot data in function of hx
    for side in sides:
        x, _, _, _, y, _, _, _, _ = data[side]
        y = y * np.power(side, (1/8))
        x = (x - 1) * side
        plt.errorbar(x, y, fmt='.', label=f'side = {side}')
    # save and show
    plt.legend(loc='upper right')
    plt.savefig(os.path.join("Plots_and_fit", title + ".png"))
    plt.show()

def plot_chi_scaling(data):
    """ Plot susceptibility scaling """

    title = "Plot scaling susceptibility Z"
    print(title + "\n")
    # axis and style
    fig = plt.figure(title)
    plt.style.use('seaborn-whitegrid')
    plt.title(title)
    plt.ylabel(r'$ \chi * side^{- \gamma / \nu} $')
    plt.xlabel('$ (hx - 1) * side  $')
    plt.xlim(-0.6, 1.5)
    plt.ylim(0.4, 2.3)
    # load and plot data in function of hx
    for side in sides:
        x, _, _, _, _, _, _, _, y = data[side]
        y = y * np.power(side, -(7/4))
        x = (x - 1) * side
        plt.errorbar(x, y, fmt='.', label=f'side = {side}')
    # save and show
    plt.legend(loc='upper right')
    plt.savefig(os.path.join("Plots_and_fit", title + ".png"))
    plt.show()

#-------------------------------------------------------------------------------

if __name__ == '__main__':

    data = load_data()

    plot_energy_gs(data)
    plot_energy_gap1(data)
    plot_energy_gap2(data)

    alpha_hx(data)

    plot_magnetization_z(data)
    plot_magnetization_x(data)
    plot_magnetization_y(data)

    plot_mag_scaling(data)
    plot_chi_scaling(data)
