"""*****************************************************************************
*
* Plot program for the outcomes of the simulation
*
*****************************************************************************"""

import os

import numpy as np
import matplotlib.pyplot as plt



#--- Contents ------------------------------------------------------------------

def load_data(side, border):
    """ Load data produced by full diagonalization """

    data = {}
    # define data file path
    filename = f"levels_{side}_hz_{0:.6f}.dat"
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
    filename = f"levels_{side}_hz_{1:.6f}.dat"
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

    binwidths = {10: 0.03125, 11: 0.03125, 12: 0.03125}

    for side in [10, 11, 12]:

        border = side*side
        binwidth = binwidths[side]
        n_bins = np.arange(0, 3.2, binwidth)

        data = load_data(side, border)

        fig, (ax0, ax1) = plt.subplots(nrows=2, figsize=(8,12))

        ax0.set_title(r'Integrable Ising | $P(s) \sim e^{-s}$')
        ax0.hist(data[0], n_bins, ls='dashed', alpha = 0.5, lw=3)

        ax1.set_title(r'Non integrable Ising | $P(s) \sim s^\gamma e^{-s^2}$')
        ax1.hist(data[1], n_bins, ls='dashed', alpha = 0.5, lw=3)

        plt.savefig(os.path.join("Plots_and_fit", f"Level spacing stat {side}.png"))
        plt.show()
