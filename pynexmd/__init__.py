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

'''
*****************************************************
PyNEXMD: Python-based Non-adiabatic EXcited-state 
         Molecular Dynamics simulations
*****************************************************

'''
__version__ = '1.0'

import os
import sys

PYNEXMD_PATH = os.getenv('PYNEXMD_PATH')

import numpy as np
from pynexmd import lib  # low-level mathematical libs (pytorch, bml, petsc, etc)
from pynexmd import mole # molecular module
from pynexmd import es   # electronic structure solver
from pynexmd import md   # molecular dynamics driver


