# Superconducting Quantum Circuit Simulations with the Harmonic Balance Method in Julia

## Description

This repository provides a set of tools for simulating superconductive quantum circuits using the harmonic balance method in Julia.

## Files Description

Here’s a summary of the key files in the repository:

- **`simulation_sweep_results.jl`**: The main code file to set up and run simulations. This script manages parameter sweeps and collects results for analysis.
  
- **`simulation_black_box.jl`**: Contains the core calculations and algorithms used in the simulations. The "black box" is where the main computation happens.
  
- **`snail_circuit.jl`**: A function to generate and simulate the Superconducting Nonlinear Asymmetric Inductive Element (SNAIL) circuit, a key component in many quantum circuit designs.

- **`utils.jl`**: Contains utility functions used throughout the project, such as data handling, custom solvers, or visualization helpers.

- **`flux_curve.txt`**: A file that maps the values of the alpha parameter of the SNAIL to the corresponding flux value to achieve the best theoretical three-wave mixing.
