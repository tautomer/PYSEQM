import numpy
from pynexmd import es, mole

mol = mole.molecule.M(
    atom = 'Li 0 0 0; H 0 0 1.75202',  # in Angstrom
    basis = 'sto3g',
    unit = 'Bohr', 
    symmetry = True,
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
