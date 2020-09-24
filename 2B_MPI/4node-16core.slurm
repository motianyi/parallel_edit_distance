#!/bin/bash

#SBATCH --error=seqalkway-tianyim-%j.err
#SBATCH --output=seqalkway-tianyim-%j.out
#SBATCH --partition=snowy
#SBATCH --time=0:40:00
#SBATCH --nodes=4-4
#SBATCH --cpus-per-task=16
#SBATCH --ntasks-per-node=1 
#SBATCH --mem=32G #SBATCH --job-name=seqalkway 
# You need to load a compiler before openmpi.
module load gcc/8.3.0
module load openmpi/3.1.4 

mpicxx -std=c++14 -fopenmp -o seqalign-mpi-tianyim seqalign-mpi-tianyim.cpp -O3
mpirun -np 4 seqalign-mpi-tianyim < mseq-big13-example.dat