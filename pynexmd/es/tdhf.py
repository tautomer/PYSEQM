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


import numpy
import os, sys

from pynexmd import lib
from pynexmd import mole
from pynexmd import __config__

class TDA(lib.StreamObject):
    '''
    CIS/TDA
    '''
    conv_tol = getattr(__config__, 'tdhf_conv_tol', 1e-9)
    nstates = getattr(__config__, 'tdhf_nstates', 3)
    singlet = getattr(__config__, 'tdhf_singlet', True)
    lindep = getattr(__config__, 'tdhf_lindep', 1e-12)
    level_shift = getattr(__config__, 'tdhf_level_shift', 0)
    max_space = getattr(__config__, 'tdhf_max_space', 50)

    def __init__(self, mf):

        self.verbose = mf.verbose
        self.stdout = mf.stdout
        self.mol = mf.mol
        self.solvent = mf.solvent
        self.mf = mf

        self.converged = None
        self.evals = None  # eigen values
        self.evecs = None  # eigenvectors

    def gen_Lxi(self, mf):
        '''
        compute Ax or 
        [ A    B][X] 
        [-B   -A][Y]
        '''
        mol = mf.mol
        mo_coeff = mf.mo_coeff
        mo_energy= mf.mo_energy
        mo_occ   = mf.mo_occ
        nocc = mf.nocc # number of occupied orbitals
        nvir = mf.nvir # number of virtual orbitals

        #1) gen V(xi) in AO
        if self.solvent is not None:
            'add solvent'
            Vxi = None
        else:
            Vxi = None


        #2) transform Vxi to MO

        hdiag = None
        #3) add diagonal term in MO



        return Vxi, hdiag

    def init_guess(self, mf):
        '''
        generate initial guess for TDA calculations
        '''

        return None

    def gen_precond(self, hdiag):

        '''
        generate preconditioner for TDA Davidson solver
        '''

        def precond(x, e, x0):
            'diagonal preconditioner'
            A = hdiag  - (e - self.level_shift)
            A[abs(A)<1.e-8] = 1.e-8
            return x / A
        return precond

class TDHF(TDA):


    '''
    '''



