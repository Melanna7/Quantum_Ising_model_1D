"""*****************************************************************************
*
* Plot program for the outcomes of the simulation
*
*****************************************************************************"""

import os

import numpy as np
import matplotlib.pyplot as plt

binwidth = 0.1
border = 1
hz_field = 0.7

#--- Contents ------------------------------------------------------------------


def load_data():
    """ Load data produced by full diagonalization """

    data = {}
    # define data file path
    filename = f"levels_hz_{0:.6f}.dat"
    file_path = os.path.join("Data_Ising", filename)
    print("Loading " + file_path)
    # load data from each file
    if os.path.isfile(file_path):
        values = np.loadtxt(file_path, unpack='True')

    data[0] = [(values[i + 1] - values[i]) for i in range(len(values) - 1)]
    data[0] = data[0][border:-border]
    norm = np.sum(data[0]) / len(data[0])
    data[0] = data[0] / norm

    # define data file path
    filename = f"levels_hz_{hz_field:.6f}.dat"
    file_path = os.path.join("Data_Ising", filename)
    print("Loading " + file_path)
    # load data from each file
    if os.path.isfile(file_path):
        values = np.loadtxt(file_path, unpack='True')

    data[1] = [abs(values[i + 1] - values[i]) for i in range(len(values) - 1)]
    data[1] = data[1][border:-border]
    norm = np.sum(data[1]) / len(data[1])
    data[1] = data[1] / norm

    return data

#-------------------------------------------------------------------------------

if __name__ == '__main__':

    data = load_data()

    n_bins = np.arange(0, 5, binwidth)

    fig, (ax0, ax1) = plt.subplots(nrows=2, figsize=(8,12))

    ax0.set_title(r'Integrable system | $P(s) \sim s^\gamma e^{-s^2}$')
    ax0.hist(data[0], n_bins, ls='dashed', alpha = 0.5, lw=3)

    ax1.set_title(r'Non integrable one | $P(s) \sim s^\gamma e^{-s^2}$')
    ax1.hist(data[1], n_bins, ls='dashed', alpha = 0.5, lw=3)

    plt.savefig(os.path.join("Plots_and_fit", "Level spacing stat.png"))
    plt.show()
