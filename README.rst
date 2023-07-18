======================
Quantum Ising model 1D
======================

In the present repository, we describe the study of the one-dimensional quantum Ising model through numerical simulation. In particular, by implementing various diagonalization techniques, the energy of the ground state and the behavior of the energy gaps between the ground state and the first two excited levels around the quantum phase transition are obtained. Furthermore, the behavior of certain physical quantities of the system around the transition is analyzed, and the effect of finite-size scaling is examined. Finally, the integrability of the transverse field Ising model is verified by examining the statistics of the level spacing.

Repository Structure
====================

The structure of the repository is as follows:

- ``class_lapacke.h`` and ``class_Eigen.h``:

  These headers, located in the ``include`` directory, contain the definition of a spin Hamiltonian operator class based on the vector and Eigen library respectively. The first implements the ZHEEV subroutine for the exact diagonalization, while the second relies on Eigen::matriX's methods. In both cases, sparse diagonalization is implemented via the ``lambda_lanczos`` library.

- ``ising_main.cpp``:

  This program calls the simulation subroutine for all the values of the transverse field and sides, collecting energy gaps, magnetization, and susceptibility measures in the ``Data_Ising`` folder. It uses the ``class_lapacke.h`` Hamiltonian, which has been proven to be faster.

- ``level_spacing.cpp``:

  This program executes the exact diagonalization of the desired spin model and saves the Eigenvalues in the ``Data_Ising`` folder. It uses the ``class_Eigen.h`` Hamiltonian, which has been proven to be more suitable for large matrix allocation and implements a random field on the longitudinal axis to break some symmetries.

- ``ising_plot_all.py`` and ``level_histogram.py``:

  These programs utilize the data in the ``Data_Ising`` folder to plot the physical quantities of interest.

- ``Heisenberg``:

  This folder contains materials that may be used to simulate the Heisenberg model properties in the next future.

- ``Tests``:

  This directory contains easy-to-use examples for testing the classes and verifying that different diagonalization methods return the same results.

- ``Plots_and_fit``:

  All produced plots are stored in this folder.

Analysis Results
================

Here are some of the plots generated from the analysis:

- Order parameter of the phase transition:

  .. image:: https://github.com/Dario-Maglio/Quantum_Ising_model_1D/blob/72e64ca7b57afae36b26f4916c3475c2212a476a/Plots_and_fit/Plot%20magnetization%20Z.png
     :align: center

- First energy gap:

  .. image:: https://github.com/Dario-Maglio/Quantum_Ising_model_1D/blob/72e64ca7b57afae36b26f4916c3475c2212a476a/Plots_and_fit/Plot%20energy%20first%20gap.png
     :align: center

- Level spacing statistics:

  .. image:: https://github.com/Dario-Maglio/Quantum_Ising_model_1D/blob/72e64ca7b57afae36b26f4916c3475c2212a476a/Plots_and_fit/Level%20spacing%20stat%2010.png
     :align: center
     :width: 80%


Feel free to explore the repository and use the provided programs for further analysis and investigation.

License
=======

This repository is licensed under the GNU General Public License v3.0 (GPL-3.0). 

See the LICENSE file for more information.

