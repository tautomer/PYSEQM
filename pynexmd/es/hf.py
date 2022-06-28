# Copyright 2022-20xx  The PyNEXMD Developers. All Rights Reserved.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
# Author: Yu Zhang zhy@lanl.gov
#

import os, sys
from pynexmd import lib
from pynexmd import mole
from pynexmd import __config__

# For code compatibility in python-2 and python-3
if sys.version_info >= (3,):
    unicode = str

def kernel(mf, conv_tol=1e-10, conv_tol_grad=None,
           dump_chk=True, dm0=None, callback=None, conv_check=True, **kwargs):

    '''
    kernel of scf driver
    '''

def energy_tot():

    return None


def get_hcore():



    return None


def get_fock():

    return None


class HF(lib.StreamObject):

    '''
    ground state SCF base class (RHF)

    verbose: int

    chkfile: str
    conv_tol : float
    max_cycle: int
    init_guess: str

    DIIS: DIIS class
    diis: 
    diis_space: init

    diis_start_crit: float (trigger diis only when error < crit)

    diis_scf: bool

    converged : bool
    etot: float
    mo_energy: orbiatl energies
    mo_occ: occupation number
    mo_coeff: orbital coefficients

    '''
    conv_tol = getattr(__config__, 'SCF_conv_tol', 1e-9)
    conv_tol_grad = getattr(__config__, 'SCF_conv_tol_grad', None)
    max_cycle = getattr(__config__, 'SCF_max_cycle', 50)

    def __init__(self, mol):

        self.mol = mol
        self.mo_energy = None
        self.mo_occ = None
        self.mo_coeff = None
        self.nocc = None
        self.nvir = None
        self.solvent = None
        self.etot = 0.0

        print('test: creating HF object')

    def scf(self):

        return None  

    def get_hcore(self):

        return None

    def get_veff(self):

        return None

    def ao2mo(self):

        return None

    def mo2ao(self):

        return None

    def get_etot(self):
        # if scf is already finised
        return self.etot
        # else run scf


class UHF(HF):

    def __init__(self, mol):
        hf.HF.__init__(self, mol)
        self.nelec = None
        self.nalpha = None
        self.nbeta = None

    def s2(self):
        'return <s^2>'

        return None

    def get_veff(mol, dm):
        '''
        unrestricted HF matrix of alpha and beta spins

        '''

