
<div align="center">
  <a href="https://github.com/e-florez/PTMC/actions/workflows/makefile.yml">
    <img src="https://github.com/e-florez/PTMC/actions/workflows/makefile.yml/badge.svg">
  </a>
  <a href="https://codecov.io/gh/username/repo">
    <img src="https://codecov.io/gh/username/repo/branch/master/graph/badge.svg">
  </a>
  <a href="https://opensource.org/licenses/MIT">
    <img src="https://img.shields.io/badge/License-MIT-yellow.svg">
  </a>
  <img src="https://img.shields.io/badge/project_status-under_development-blue">
</div>

# Parallel Tempering Monte Carlo (PTMC)
## Introduction

Monte Carlo simulations are a powerful tool to study the thermodynamic properties of atomic and molecular clusters, including their phase transitions. This project presents a new Fortran (2003) code for simulating phase transitions using the Parallel Tempering Monte Carlo (PTMC) method.

## Description

The PTMC method is an improvement over previous methods due to its ability to efficiently explore the phase space and accurately compute the interaction energy using many-body expansion techniques such as MP2 and coupled-cluster theory for atomic and molecular clusters. We apply the PTMC method to study the melting transitions of noble gases (Ne, Ar, Kr, and Xe) and water clusters, as well as a novel study of melting transitions of neon clusters in a homogeneous ultra-high magnetic field of up to 70,500 Tesla.

## Requirements

- GNU Make 4.3
- gcc 11.3.0
- Open MPI 4.1.2

## Installation

1. Clone the repository from GitHub:

    `git clone https://github.com/e-florez/PTMC.git`

2. Change to the repository directory:

    `cd PTMC`

3. Run the `configure` script to generate the `Makefile`:

    `./configure --src-dir=src`

4. Compile the code using the generated `Makefile` 

    `make [option]`

    with one of the following options:

- `make`: Compile the code with sequential execution (single CPU).
- `make openmp`: Compile the code with OpenMP parallelization (multi-threading).
- `make array`: Compile the code with array bounds checking enabled (single CPU).
- `make debug`: Compile the code with debugging options enabled (single CPU).
- `make debug-openmp`: Compile the code with debugging options enabled and OpenMP parallelization (multi-threading).


5. To clean up the generated files and executables, run:

    `make clean`


## Usage

Run the compiled executable with the desired options:

`./PTMC.exe [.inp file]`


## License

This project is licensed under the MIT License.

## Project status

This project is currently under development.
