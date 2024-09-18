from ase.io import read, write 
from ase import atoms 
from ase.io.vasp import read_vasp
from ase.visualize import view
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
import random
import numpy as np
from ase.md import MDLogger
from ase.md.npt import NPT
from ase import units
from mace.calculators import mace_mp
from tqdm import tqdm
from ase.md.analysis import DiffusionCoefficient
import os
import matplotlib.pyplot as plt

def get_index_arr(trajectory, atom='Li', start=0):
    # takes in a trajectory and identifies the indices of the given atom in frame start
    #
    # trajectory - ASE trajectory, list of ASE atoms objects
    # atom - string element symbol of element you wish to examine
    # start - integer, frame to examine - useful if you wish to skip the first few steps of the trajectory 
    #
    # returns
    # index_arr - list, list of all indices, for ASE atoms object, of a given element in a given trajectory
    index_arr=[]
    for i, j in enumerate(trajectory[start].get_chemical_symbols()):
        if j == atom:
            index_arr.append(i)
    return index_arr


def get_diffusion_coeff(trajectory, timestep, atom='Li', start=0):
    # lanches the ASE diffusion coefficient method
    #
    # trajectory - ASE trajectory, list of ASE atoms objects
    # timestep - float, the timestep of the molecular dynamics trajectory
    # atom - string element symbol of element you wish to examine
    # start - integer, frame to examine - useful if you wish to skip the first few steps of the trajectory
    #
    # returns
    # returns the difusion coefficient in cm^2/S
    index_arr=get_index_arr(trajectory, atom, start)

    diff = DiffusionCoefficient(trajectory, timestep=timestep, atom_indices=index_arr)

    coeff = diff.get_diffusion_coefficients()
    return coeff[0][0] * (10**-1) * units.fs


def run_md(temp_arr=[525], repeats=5, ts=2*units.fs, tts=80*units.fs, step=100, model='medium', floats='float32'):

    # given a set of inputs, sets up and runs a range of MACE-MP-0 MD simmulations !!!! in series !!!!
    # temp_arr - list of integers, a list of temperatures you wish to run MD for
    # repeats - number of MD simulations that should be run at each temperature
    # ts - float, timestep usually in fs
    # tts - characteristic time scale for the nose hoover thermostat usually in fs. A good starting point is 2*ts
    # step - integer, number of steps for each MD run
    # model - string, controls the MACE-MP-0 model that will be used, choice of small, medium and large
    # floats - string, controls what float type the MACE-MP-0 calculator will use, float32 or float64
    for temp in temp_arr:
        for run in range(repeats):

            ## set random seed - new for each run repitition
            random.seed(run)

            ## read in init structure
            init_conf = read('g_li3ps4.xyz')

            ## correct rounding error
            a = init_conf.get_cell()
            a[1,0]=0.0;a[2,0]=0.0;a[2,1]=0.0
            init_conf.set_cell(a)

            ## launch and set mace calculator
            calc = mace_mp(model=model, device='cpu', default_dtype=floats)
            init_conf.calc = calc

            ## define write file format
            def write_frame():
                dyn.atoms.write('mace_md_'+str(temp)+'K_repeat_'+str(run)+'.xyz', append=True)

            ## initialise the MD
            MaxwellBoltzmannDistribution(init_conf, temperature_K=temp)
            dyn = NPT(init_conf, externalstress=1., timestep =ts*units.fs, ttime=tts*units.fs, temperature_K=temp)
            dyn.attach(write_frame, interval=1)

            ## run MD
            dyn.run(step)



def extract_diff_coeff(temp_arr=[525], repeats=5, atom='Li', start=0, ts=2*units.fs):
    # reads in the final results and extracts the diffusion coefficient
    # 
    # temp_arr - list of integers, a list of temperatures you wish to run MD for
    # repeats - number of MD simulations that should be run at each temperature
    # atom - string element symbol of element you wish to examine
    # start - integer, frame to examine - useful if you wish to skip the first few steps of the trajectory
    # ts - float, timestep usually in fs
    # !!!! make sure these inputs as the same as used when running the MD !!!!
    #
    # returns 
    #dict_diff_coeff - a dictionary with each temperature as a key and the diffusion constants at each repeat in an array
    dict_diff_coeff = {}
    for temp in temp_arr:
        dict_diff_coeff[str(temp)]=[]
        for run in range(repeats):
            trajectory = read('mace_md_'+str(temp)+'K_repeat_'+str(run)+'.xyz', index = ':')
            coeff = get_diffusion_coeff(trajectory, ts, atom='Li', start=0)
            dict_diff_coeff[str(temp)].append(coeff)
            
    return dict_diff_coeff


def plot_diff_coeff(dict_diff_coeff):
    # given a dictionary of diffusion coefficients, plot all values using a log y scale and change the x scale to be 1000/temperature
    #
    # dict_diff_coeff - dictionary of arrays of floats. Keys correspond to examined temperatures and values correspond to diffusion constants at repeats
    #
    # return
    # coefficients1 - float y = mx + c values
    # cov - 2x2 matrix of floats, covariance matrix
    # also returns the plot for Lithium diffusion 
    temp_arr2=[]
    y=[]
    for i in dict_diff_coeff:
        for j in dict_diff_coeff[i]:
            temp_arr2.append(int(i))
            y.append(float(j))
    x = []
    for i in temp_arr2:
        x.append(1000/i)
        
        
    coefficients1,cov = np.polyfit(x, y, 1, cov=True)
    print(coefficients1,cov)
    poly2 = np.poly1d(coefficients1)
    best_fit_line3 = poly2(x)
    
    best_fit = []
    for i in x:
        best_fit.append(coefficients1[0]*i + coefficients1[1])
    # axs[1, 0].plot(Tinv_small, best_fit_line3, label='vasp try',color = 'blue', linestyle='-')
    plt.plot(x, best_fit,color = 'blue', linestyle='dotted')    
        
        
        
        
    plt.scatter(x, y)
#     plt.yscale('log')
    plt.xlabel('1000 * 1/T K^-1', fontsize=24)
    plt.ylabel('D* [cm2/S]', fontsize=24)
    plt.title('Lithium diffusion in g-Li3PS4')
    plt.show()
    return coefficients1,cov

def plot_average_diff_coeff(dict_diff_coeff):
    # given a dictionary of diffusion coefficients, plot average diffusion constant at a temperature using a log y scale and change the x scale to be 1000/temperature
    #
    # dict_diff_coeff - dictionary of arrays of floats. Keys correspond to examined temperatures and values correspond to diffusion constants at repeats
    #
    # return
    # coefficients1 - float y = mx + c values
    # cov - 2x2 matrix of floats, covariance matrix
    # also returns the plot for Lithium diffusion 
    y=[] # average diffusion constant for a temperature
    y_full=[] # all diffusion constants, including repeats
    x = [] # temperatures examined, len(x) = len(y)
    x_full=[] # temperatures examined but 1 entry for each repeat at each entry, len(x_full) = len(y_full)
    y_mat =[] # contains all diffusions constants but each temperature has its own entry

    # make 4 arrays, x arrays containt temperatures consistent to y and y_full values respectively
    for i in dict_diff_coeff:
        x.append(1000/int(i))
        for j in range(len(dict_diff_coeff[i])):
            x_full.append(1000/int(i))
            y_full.append(dict_diff_coeff[i][j])
        y.append(np.average(dict_diff_coeff[i]))
        y_mat.append(dict_diff_coeff[i])

    coefficients1,cov = np.polyfit(x_full, y_full, 1, cov=True)
    poly2 = np.poly1d(coefficients1)
    best_fit_line3 = poly2(x)
    plt.plot(x, best_fit_line3,color = 'blue', linestyle='dotted')

# Calculate mean and standard deviation
    y_std_devs =[]
    for i in y_mat:
        y_std_devs.append(np.std(i))

# Plot
    plt.errorbar(x, y, yerr=y_std_devs, fmt='o', capsize=5)
    plt.title('Mean Values with Error Bars at Different Temperatures')






    plt.scatter(x, y)
    plt.xlabel('1000 * 1/T K^-1', fontsize=24)
    plt.ylabel('D* [cm2/S]', fontsize=24)
    plt.title('Lithium diffusion in g-Li3PS4')
    plt.show()
    return coefficients1,cov


def room_temp_diff_coeff(coefficients, cov, temp=300):
    # extrapolates the line of bestfit to find diffusion at a given temperature
    # with errors
    # coefficients - float y = mx + c values
    # cov - 2x2 matrix of floats, covariance matrix
    # temp - integer, the temperature at which you wish to estimate the diffusion constant at
    #
    # returns
    # y - float, the diffusion constant at a given temperature
    # y_err - float, the error in the diffusion constant at given temperature

    inv_temp = 1000/temp

    m,c = coefficients
    m_err, c_err = np.sqrt(np.diag(cov))

    y = (m* inv_temp) + c
    y_err = np.sqrt(m_err**2 + (inv_temp * c_err)**2)
    return y, y_err

def run_md_bash(structure_file, temp_arr=[525], repeats=5, ts=2*units.fs, tts=80*units.fs, step=100, model='medium', floats='float32', device='cpu'):
    # set up and run the MD !!!! in paralel !!!! using a bash script
    # structure - string, filename of the starting structure for the MD runs. in a format that is supported by ASE read(). '.xyz' is recommended
    # temp_arr - list of integers, a list of temperatures you wish to run MD for
    # repeats - number of MD simulations that should be run at each temperature
    # ts - float, timestep usually in fs
    # tts - characteristic time scale for the nose hoover thermostat usually in fs. A good starting point is 2*ts
    # step - integer, number of steps for each MD run
    # model - string, controls the MACE-MP-0 model that will be used, choice of small, medium and large
    # floats - string, controls what float type the MACE-MP-0 calculator will use, float32 or float64
    # device - string, controls what device the calculations should be run on, cpu or gpu

    f = open("run_md.sh", "w")
    f.write("#!/bin/sh\n")
    for temp in temp_arr:
        for run in range(repeats):
            f.write('python mace_nvt.py '+str(structure_file)+' '+str(temp)+' '+str(run)+' '+str(ts)+' '+str(tts)+ ' '+str(step)+' '+str(model)+' '+str(floats)+' '+str(device)+' &\n')
    f.write('wait')
    f.close()
    os.system('./run_md.sh >out.txt 2> out.txt')

def diffusion_plot(structure_file,temp_arr=[525], repeats=5, ts=2*units.fs, tts=80*units.fs, step=100, atom='Li', start=0, model='medium', floats='float32', device='cpu'):
    # given a set of inputs, sets up and runs a range of MACE-MP-0 MD simmulations, extract the coefficients and plot the average
    #
    # structure - string, filename of the starting structure for the MD runs. in a format that is supported by ASE read(). '.xyz' is recommended
    # temp_arr - list of integers, a list of temperatures you wish to run MD for
    # repeats - number of MD simulations that should be run at each temperature
    # ts - float, timestep usually in fs
    # tts - characteristic time scale for the nose hoover thermostat usually in fs. A good starting point is 2*ts
    # step - integer, number of steps for each MD run
    # model - string, controls the MACE-MP-0 model that will be used, choice of small, medium and large
    # floats - string, controls what float type the MACE-MP-0 calculator will use, float32 or float64
    # device - string, controls what device the calculations should be run on, cpu or gpu
    #
    # returns
    # coefficients1 - float y = mx + c values
    # cov - 2x2 matrix of floats, covariance matrix
    run_md_bash(structure_file, temp_arr, repeats, ts, tts, step, model, floats, device)
    dict_diff_coeff = extract_diff_coeff(temp_arr, repeats, atom, start, ts)
    coefficients,cov = plot_average_diff_coeff(dict_diff_coeff)
    return coefficients, dict_diff_coeff,cov
