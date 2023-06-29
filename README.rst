Quantum Ising model 1D
======================

In the present repository, we study the one-dimensional Quantum Ising model
at 0°K via numerical simulation. In particular, we analyze the behaviour of the
system around the phase transition and verify the effect of finite-size scaling
effect. In doing so, we take advantage of the Eigen library diagonalization
procedures and the beta-Lanczos repository for dealing with sparse Hamiltonians.

The scheme of the repository is the following:

- class_Eigen.h contains the Hamiltonian class and the HamiltParameters structure. The class is based on the Eigen library and construct a dense or sparse Hamiltonian depending on the value of a sparse_flag. It also provides a diagonalization method from the Eigen library in the dense case and a Lanczos algorithm for the sparse one.

- class_lapacke.h contains another implementation of the Hamiltonian class based on the vector library. It provides a diagonalization method via LAPACKE in the dense matrix case.

- ising_main.cpp and heise_main.cpp uses the class_Eigen.h file to instantiate different Hamiltonians varying the side and the fields values. It also collects data for the energy gaps and magnetization values in the Data_* folder.

- *_plot_all.py produces the plot in the Plots_* folders.

Easy to use examples are given in the test folder.
