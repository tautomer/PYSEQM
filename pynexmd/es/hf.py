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

    '''

class UHF(HF):


    return None


