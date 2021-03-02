# Two-qubit CZ gate implementation with trapped neutral atoms: a numerical simulation

## Abstract

Multi-qubit gates can be implemented with trapped neutral atoms by driving them to highly excited Rydberg states, in which nearby Rydberg atoms can interact very strongly and in a controlled way. In this Report, we simulate the dynamics of a system under the application of a controlled-phase gate physically implemented in [[1](https://arxiv.org/abs/1908.06101)]. After a physical and theoretical description of the two-qubit CZ gate design, the dynamics of a two- qubit system is simulated in both the perfect and imperfect Rydberg blockade regime by varying the optimal parameters suggested. Then, we consider the CZ gate in the case of a chain of N-qubit. A Gaussian noise is introduced on the optimal parameters to investigate noise effects with an increasing number of qubits. Numerical methods for performing unitary time evolution are also discussed and their performances are tested.

## Contents

1. Introduction
2. Physical implementation of single-qubit and two-qubit gates
3. Theoretical design of two-qubit CZ gate
4. Code implementation of two-qubit CZ gate
    1. Two-qubit system
    2. Chain of N-qubit
5. Numerical methods for time-dependent Schrodinger equation solvers
    1. Spectral method
    2. Crank-Nicolson method
    3. Crank-Nicolson method with LU decomposition
6. Results
    1. CZ gate in a two-qubit system
        1. Correct behaviour check and optimal parameters variation
        2. Measurement process with Gaussian noise on τ and ∆
    2. CZ gate in a chain of N-qubit
    3. Timing analysis of time evolution methods
7. Conclusions


## Authors

* [**Alice Pagano**](https://github.com/AlicePagano) (University of Padua)
* [**Michele Puppin**](https://github.com/michelepuppin) (University of Padua)

## Content of the folder

The repository is organized as follows:
* **`code`**: folder with all the source code of the project, in particular:
    * **`physical_behavior`**: jupyter notebooks to implement the CZ gate, test its correctness and investigate noise effects;
    * **`timing_analysis`**: jupyter notebooks to test and plot timing analysis results;
* **`report`**: folder with a report of the project, including also the `.tex` source files.
