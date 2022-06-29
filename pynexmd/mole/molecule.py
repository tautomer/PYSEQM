
import os, sys
from os import path
import torch
import logging
import numpy
from pynexmd.data import elements
from pynexmd.data.elements import ELEMENTS, ELEMENTS_PROTON, \
        _rm_digit, atom_charge, _symbol, _std_symbol, _atom_symbol, is_ghost_atom, \
        _std_symbol_without_ghost
from pynexmd.lib import StreamObject

from pynexmd.seqm.seqm_functions.constants import Constants, ev_kcalpmol
from pynexmd.seqm.seqm_basics import  Parser, Hamiltonian, Pack_Parameters, Energy

# For code compatibility in python-2 and python-3
if sys.version_info >= (3,):
    unicode = str

# copy from pyscf
def M(**kwargs):
    r'''This is a shortcut to build up Mole object.

    Args: Same to :func:`Mole.build`

    Examples:

    >>> from pynexmd import mole
    >>> mol = mole.M(atom='H 0 0 0; F 0 0 1', basis='6-31g')
    '''
    mol = molecule()
    mol.build(**kwargs)
    return mol


def format_atom(atoms, origin=0, axes=None,
                unit='Ang'):
    '''Convert the input :attr:`Mole.atom` to the internal data format.
    Including, changing the nuclear charge to atom symbol, converting the
    coordinates to AU, rotate and shift the molecule.
    If the :attr:`~Mole.atom` is a string, it takes ";" and "\\n"
    for the mark to separate atoms;  "," and arbitrary length of blank space
    to spearate the individual terms for an atom.  Blank line will be ignored.

    Args:
        atoms : list or str
            the same to :attr:`Mole.atom`

    Kwargs:
        origin : ndarray
            new axis origin.
        axes : ndarray
            (new_x, new_y, new_z), new coordinates
        unit : str or number
            If unit is one of strings (B, b, Bohr, bohr, AU, au), the coordinates
            of the input atoms are the atomic unit;  If unit is one of strings
            (A, a, Angstrom, angstrom, Ang, ang), the coordinates are in the
            unit of angstrom;  If a number is given, the number are considered
            as the Bohr value (in angstrom), which should be around 0.53.
            Set unit=1 if wishing to preserve the unit of the coordinates.

    Returns:
        "atoms" in the internal format. The internal format is
            | atom = [[atom1, (x, y, z)],
            |         [atom2, (x, y, z)],
            |         ...
            |         [atomN, (x, y, z)]]

    Examples:

    >>> gto.format_atom('9,0,0,0; h@1 0 0 1', origin=(1,1,1))
    [['F', [-1.0, -1.0, -1.0]], ['H@1', [-1.0, -1.0, 0.0]]]
    >>> gto.format_atom(['9,0,0,0', (1, (0, 0, 1))], origin=(1,1,1))
    [['F', [-1.0, -1.0, -1.0]], ['H', [-1, -1, 0]]]
    '''
    def str2atm(line):
        dat = line.split()
        try:
            coords = [float(x) for x in dat[1:4]]
        except ValueError:
            if DISABLE_EVAL:
                raise ValueError('Failed to parse geometry %s' % line)
            else:
                coords = list(eval(','.join(dat[1:4])))
        if len(coords) != 3:
            raise ValueError('Coordinates error in %s' % line)
        return [_atom_symbol(dat[0]), coords]

    if isinstance(atoms, (str, unicode)):
        # The input atoms points to a geometry file
        if os.path.isfile(atoms):
            try:
                atoms = fromfile(atoms)
            except ValueError:
                sys.stderr.write('\nFailed to parse geometry file  %s\n\n' % atoms)
                raise

        atoms = str(atoms.replace(';','\n').replace(',',' ').replace('\t',' '))
        fmt_atoms = []
        for dat in atoms.split('\n'):
            dat = dat.strip()
            if dat and dat[0] != '#':
                fmt_atoms.append(dat)

        if len(fmt_atoms[0].split()) < 4:
            fmt_atoms = from_zmatrix('\n'.join(fmt_atoms))
        else:
            fmt_atoms = [str2atm(line) for line in fmt_atoms]
    else:
        fmt_atoms = []
        for atom in atoms:
            if isinstance(atom, (str, unicode)):
                if atom.lstrip()[0] != '#':
                    fmt_atoms.append(str2atm(atom.replace(',',' ')))
            else:
                if isinstance(atom[1], (int, float)):
                    fmt_atoms.append([_atom_symbol(atom[0]), atom[1:4]])
                else:
                    fmt_atoms.append([_atom_symbol(atom[0]), atom[1]])

    if len(fmt_atoms) == 0:
        return []

    if axes is None:
        axes = numpy.eye(3)

    if isinstance(unit, (str, unicode)):
        if unit.upper().startswith(('B', 'AU')):
            unit = 1.
        else: #unit[:3].upper() == 'ANG':
            unit = 1./param.BOHR
    else:
        unit = 1./unit

    c = numpy.array([a[1] for a in fmt_atoms], dtype=numpy.double)
    c = numpy.einsum('ix,kx->ki', axes * unit, c - origin)
    z = [a[0] for a in fmt_atoms]
    return list(zip(z, c.tolist())), c.tolist()


class molecule(StreamObject):
    '''
    the basic class to hold 
    1) molecular structure
    2) gloabl options

    Attributes:
       verbose: init

       charge : int
       spin   : int or None
       atom : list
              format of atom:
       unit : str: au [default]
       output : str: 
       basis : str ['AM1', 'PM6']

    '''

    def __init__(self, **kwargs):

        self.output = None

        self.verbose = 0
        self.charge = 0
        self.spin = 0 # 2j = n_alpha - n_beta
        self.symmetry = False
        self.basis = 'sto-3g' # not used in semiempirical

        self.atom = [] 

        self.basis = None
        self.theory= 'AM1' # [DFT(NWchem), DFTB]
        self.unit = 'Ang'
        self.symmetry = False
        self.device = None

    def build(self, atom=None,
            basis = None,
            theory=None, #'AM1', DFTB, 'DFT' (using pySCF),
            unit=None,
            symmetry=None,
            verbose = None,
            charge=0, spin=0,
            device = None):


        if verbose is not None: self.verbose = verbose
        if unit is not None: self.unit = unit
        if theory is not None: self.theory = theory
        if basis is not None: self.basis = basis
        if atom is not None: self.atom = atom
        if charge is not None: self.charge = charge
        if spin != 0: self.spin = spin
        if device is not None: self.device = device

        self._atom, coordinates = self.format_atom(self.atom, unit=self.unit)
        print('test-atom=', self._atom)

        self.nats = len(self._atom)
        self.atom_order = [atom_charge(atom[0]) for atom in self._atom]
        
        # TODO, charge for AM1 is different (as only valence electrons are considered)
        self.nelectrons = sum(atom_charge(atom[0]) for atom in self._atom)

        self.const = Constants().to(device)
        self.species = torch.as_tensor([self.atom_order], dtype=torch.int64, device=self.device)
        if self.theory == 'AM1':
            tore = self.const.tore
            self.nelectrons = torch.sum(tore[self.species],dim=1).reshape(-1).type(torch.int64)
        self.nelectrons -= self.charge

        self.coordinates = torch.tensor([coordinates], device=self.device)
        self.coordinates.requires_grad_(True)

        print('test-total electrons=', self.nelectrons)
        print('test-number of atoms', self.nats)
        print('test-atom_order=', self.atom_order)
        print('test-species=', self.species)
        print('test-coordinates=', coordinates)
        
        here = path.abspath(path.dirname(__file__))
        print('here=', here)

        self.seqm_parameters = {
                   'method' : 'AM1',  # AM1, MNDO, PM#
                   'scf_eps' : 1.0e-6,  # unit eV, change of electric energy, as nuclear energy doesnt' change during SCF
                   'scf_converger' : [2,0.0], # converger used for scf loop
                                         # [0, 0.1], [0, alpha] constant mixing, P = alpha*P + (1.0-alpha)*Pnew
                                         # [1], adaptive mixing
                                         # [2], adaptive mixing, then pulay
                   'sp2' : [True, 1.0e-5],  # whether to use sp2 algorithm in scf loop,
                                            #[True, eps] or [False], eps for SP2 conve criteria
                   'elements' : elements, #[0,1,6,8],
                   'learned' : ['U_ss'], # learned parameters name list, e.g ['U_ss']
                   'parameter_file_dir' : here+'/../lib/params/', # file directory for other required parameters
                   'pair_outer_cutoff' : 1.0e10, # consistent with the unit on coordinates
                }

        #parser is not needed here, just use it to get Z and create "learned" parameters
        #prepare a fake learned parameters: learnedpar
        self.parser = Parser(self.seqm_parameters).to(device)
        self.nmol, self.molsize, \
        self.nHeavy, self.nHydro, self.nocc, \
        self.Z, self.maskd, self.atom_molid, \
        self.mask, self.pair_molid, self.ni, self.nj, \
        self.idxi, self.idxj, self.xij, self.rij = self.parser(self.const, self)

        print('test-nocc', self.nocc)


        return self

    kernel = build

    def format_atom(self, atom, origin=0, axes=None, unit='Ang'):
        return format_atom(atom, origin, axes, unit)

    def gen_int1e(self):

        if self.device == 'torch':
            return 
        elif self.device == 'BML':
            return None
        else:
            'using numpy'
            return None


