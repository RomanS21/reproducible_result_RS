# reproducible_result_RS
This program is designed to run MACE based NVT MD and extract the diffusion coefficients, when done for multiple temperatures, this is then plotted. From this an estimate of the diffusion constant can be calculted with proppogated errors. This program takes in a structure file like a .xyz, and applies a MACE model. The model can be a MACE-MP-0 model or a custom potential, there are also extra parameters that are carried over from how you would launch MACE calculations from ASE. You can then specify an array of temperatures, the number of repeats and parameters for the NVT dynamics carried over from the ASE Nose-Hoover implementation. The program then runs the calculations, in series or in a crude parallel verison, on cpu or gpu and saves the output trajectories. These can then be read in and the diffusion constants extracted, an option is given to ignore a certain number of initial structures for poorly converged MD runs. There are also scripts to plot the values, either all of them or averaged for each temperature, from this a line of best fit is calculated and a function has been written to use this to extrapolate the diffusion constant at a given temperature with propogated errors in this value. 

In this mini tutorial, you will be intoruced to MACE-MP-0 a foundation machine learned interatomic potential. With this potential, you will run NVT molecular dynamics for g-Li3PS4 at various temperatues, with multiple repeats. The trajectories will then be analysed to calculated a mean sqaure displacement (MSD) the gradient of which corresponds to the difussion constant. The diffusion constant is the value we are interested in. Having obtained multiple diffusion constants for given temperatures, a mean value will be calculated and plotted at each temperature. These values will then be used to extrapolate the diffusion constant at room temperature. Since we have multiple repeats at each temperature we can approximate some of the error in our calculated room temperature diffusion constant and compare it to a precalculated value.

Please proceed with installing the following dependencies that have been made based on the assumption that hetmathsys is being used. Once you are done, please launch the jupyter notebook. Finally, most functions can be found in the utils folder, excepting the mace_nvt.py script that runs the MACE MD itself.

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


