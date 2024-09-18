import sys
import os
from ase.io import read, write 
from ase import atoms 
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
import random
import numpy as np
from ase.md import MDLogger
from ase.md.npt import NPT
from ase import units
from mace.calculators import mace_mp, MACECalculator
from ase.md.analysis import DiffusionCoefficient



## read in the command line inputs from the bash run script
init_file,temp,run,ts,tts,step,model,floats,device = sys.argv[1:]

## set the random seed unique for each repeat
random.seed(run)

print('initialise')            
## read in init structure
init_conf = read(init_file)
            
## When reading the VASP POSCAR file, later saved as .xyz, 0.0 can be read as an infenitesimal number correct this to be 0 as NVT dynamics requires upper triangle
a = init_conf.get_cell()
a[1,0]=0.0;a[2,0]=0.0;a[2,1]=0.0
init_conf.set_cell(a)
            
## launch and set mace calculator, if small, medium or large is specified then the use of mace mp is implied, else a path to the mace model is read in a used
mp_models = ['small','medium','large']
if model.strip().lower() in mp_models:
    calc = mace_mp(model=model, device=device, default_dtype=floats)
else:
    calc = MACECalculator(model_paths=[model], device=device, default_dtype=floats)
init_conf.calc = calc

## remove previous repeats, this avoids any append to structures, this can be suppressed but for the purpose of this result it could make the process more confusing
os.system('rm mace_md_'+str(temp)+'K_repeat_'+str(run)+'.xyz')

## define write file format 
def write_frame():
    dyn.atoms.write('mace_md_'+str(temp)+'K_repeat_'+str(run)+'.xyz', append=True)
            
## initialise the MD
MaxwellBoltzmannDistribution(init_conf, temperature_K=float(temp))
dyn = NPT(init_conf, externalstress=1., timestep =float(ts), ttime=float(tts), temperature_K=float(temp))
dyn.attach(write_frame, interval=1)
            
## run MD
dyn.run(int(step))

