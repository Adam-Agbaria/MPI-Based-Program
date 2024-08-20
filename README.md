# MPI Point Averaging with Cartesian Topology

This project implements a parallel algorithm using MPI (Message Passing Interface) to distribute and process points in a 2D space. The goal is to iteratively adjust the positions of these points based on their neighbors until all points converge within a specified distance `D` from the origin or until a maximum number of iterations is reached. The program uses a Cartesian topology to structure the MPI processes.

## Table of Contents
- [Features](#features)
- [Requirements](#requirements)
- [Installation](#installation)
- [Usage](#usage)


## Features

- **Parallel Point Distribution**: Distributes points across MPI processes.
- **Cartesian Topology**: Organizes processes in a 2D Cartesian grid for neighbor communication.
- **Iterative Point Adjustment**: Iteratively adjusts each point's coordinates based on its neighbors.
- **Convergence Check**: Terminates if all points converge within a specified distance `D`.
- **Flexible Input**: Reads points and configuration parameters from a text file.
- **Gather Results**: Collects and outputs the final coordinates of all points.

## Requirements

- **MPI**: MPI library (e.g., OpenMPI, MPICH).
- **C Compiler**: GCC or another C compiler that supports MPI.
- **Make**: (Optional) For building the project.

## Build the Project

1. **Compile the program** using `mpicc`:
   ```bash
   mpicc -o mpi_point_avg main.c -lm
## Prepare Input File

Ensure you have an input file named `data.txt` in the same directory as the executable (see [Input File Format](#input-file-format) for details).

## Usage

Run the program using `mpirun` or `mpiexec`, specifying the number of processes:
  ```bash
  mpirun -np <number-of-processes> ./mpi_point_avg

Example:
mpirun -np 4 ./mpi_point_avg
This will run the program with 4 MPI processes.
