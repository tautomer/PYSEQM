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
            mol,      # molecule 
            esdriver, # electronic structure driver
            *args,
            ):
        """
        esdriver: chose the driver for electronic structure:
                semi-empirical, DFT, CCSD (pySCF), QED-CCSD, HIPNN, etc.
        args: dicitonary of other arguments: 
              1) classical time step
              2) thermostat
              3) number of classical steps
              ...
        """


        self.mol = mol
        self.esdriver = esdriver
        #....

    def verlet(self):

        return None

    def thermostat(self, vel):
        '''
        thermostate, rescale velocities
        '''

        return None

    def velocity_rescale(self, vec, deltaE):
        """
        energy difference before and after rescaling
        vec: direciton of velocity rescaling
        """
        return None

    def output(self):

        return None
