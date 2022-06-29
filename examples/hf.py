import os, sys
import numpy
import torch
from pynexmd import es, mole

#import logging
#logging.basicConfig(level=logging.INFO)

torch.set_default_dtype(torch.float64)

if torch.cuda.is_available():
    device = torch.device('cuda')
    print('GPU is used')
else:
    device = torch.device('cpu')
    print('CPU/numpy is used')



mol = mole.molecule.M(
    atom = 'Li 0 0 0; H 0 0 1.75202',  # in Angstrom
    unit = 'Bohr',
    theory = 'AM1',
    device = device
)


hf = es.hf.HF(mol)
hf.max_cycle = 200
hf.conv_tol = 1.e-8
hf.diis_space = 10

'''
#hf.polariton = True
#hf.gcoup = 0.01
mf = hf.run()
print('hf energy',mf.hf_energy)
'''
