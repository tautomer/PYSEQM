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

#import torch
import os, sys
from pynexmd import lib
from pynexmd import mole
from pynexmd import __config__


class MD(lib.StreamObject):
    '''
    '''
    def __init__(self, 
            theory, 
            coordinates, 
            seqm_paramters = None,
            init_vec = None,
            *args,
            ):
        """
        theory: chose the theory level for electronic structure:
                semi-empirical, DFT, CCSD (pySCF)
        coordinates: coordinates of the system of interest
        mdtype: "SH" (default), "AIMC", "Ehrenfest"
        args: dicitonary of other arguments: 
              1) number of excited states
              2) number of classical steps
              3) mdtypes 
              4) classical time step
              5) classic/quantum time ration
              6) decoherence type
              7) crossing check
              ...
              20) verbose
        """


        self.theory = theory
        self.coords = coordinates
        self.init_vec = init_vec
        
        #self.mdtype = mdtype

    def verlet_nve(self):

        return None

    def verlet_nvt(self):

        return None

    def velocity_rescale(self, vec, deltaE):
        """
        energy difference before and after rescaling
        vec: direciton of velocity rescaling
        """
        return None

    def output(self):

        return None
