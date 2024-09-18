# reproducible_result_RS
Reproducible result for Roman Shantsila, PX915
# Dependancy 
For this workshop you will need access to Python, PyTorch, MACE and ASE

Python  >= 3.7

PyTorch >= 1.12



Python:

Python 3.10.4 was used on HetMathsys for the calculations to load this use

ml GCCcore/11.3.0 Python/3.10.4

MACE:

MACE can be installed as directed here https://github.com/ACEsuit/mace

Python  >= 3.7

PyTorch >= 1.12

Install PyTorch here https://pytorch.org/get-started/locally/

For pip install of MACE run

pip install --upgrade pip

pip install mace-torch

ASE:

ASE 3.22.1 is recomended 

For pip install run

pip install ASE==3.22.1

This tutotial will also use Python modules numpy, matplotlib, os, and random. These are assumed to be included with your loaded Python.


In this mini tutorial, you will be intoruced to MACE-MP-0 a foundation machine learned interatomic potential. With this potential, you will run NVT molecular dynamics for g-Li3PS4 at various temperatues, with multiple repeats. The trajectories will then be analysed to calculated a mean sqaure displacement (MSD) the gradient of which corresponds to the difussion constant. The diffusion constant is the value we are interested in. Having obtaiend multiple diffusion constants for given temperatures, a mean value will be calculated and plotted at each temperature. These values will then be used to extrapolate the diffusion constant at room temperature. Since we have multiple repeats at each temperature we can approximate some of the error in our calculated room temperature diffusion constant and compare it to a precalculated value.
