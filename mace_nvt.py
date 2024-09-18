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




init_file,temp,run,ts,tts,step,model,floats,device = sys.argv[1:]
random.seed(run)
print('initialise')            
## read in init structure
init_conf = read(init_file)
            
## correct rounding error
a = init_conf.get_cell()
a[1,0]=0.0;a[2,0]=0.0;a[2,1]=0.0
init_conf.set_cell(a)
            
## launch and set mace calculator
mp_models = ['small','medium','large']
if model.strip().lower() in mp_models:
    calc = mace_mp(model=model, device=device, default_dtype=floats)
else:
    calc = MACECalculator(model_paths=[model], device=device, default_dtype=floats)
init_conf.calc = calc
            
os.system('rm mace_md_'+str(temp)+'K_repeat_'+str(run)+'.xyz')
## define write file format 
def write_frame():
    dyn.atoms.write('mace_md_'+str(temp)+'K_repeat_'+str(run)+'.xyz', append=True)
            
## initialise the MD
MaxwellBoltzmannDistribution(init_conf, temperature_K=float(temp))
dyn = NPT(init_conf, externalstress=1., timestep =float(ts)*units.fs, ttime=float(tts)*units.fs, temperature_K=float(temp))
dyn.attach(write_frame, interval=1)
            
## run MD
dyn.run(int(step))

