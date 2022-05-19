#

"""

"""

import os, sys
import numpy as np
import torch

#
# the idea from pyscf
class StreamOjbect(object):

    verbose = 0
    stdout = sys.stdout

    def kernel(self, *args, **kwargs):
        '''
        Kernel function is the main driver of a method.  Every method should
        define the kernel function as the entry of the calculation.  Note the
        return value of kernel function is not strictly defined.  It can be
        anything related to the method (such as the energy, the wave-function,
        the DFT mesh grids etc.).
        '''
        pass

    def pre_kernel(self, envs):
        '''
        A hook to be run before the main body of kernel function is executed.
        Internal variables are exposed to pre_kernel through the "envs"
        dictionary.  Return value of pre_kernel function is not required.
        '''
        pass

    def post_kernel(self, envs):
        '''
        A hook to be run after the main body of the kernel function.  Internal
        variables are exposed to post_kernel through the "envs" dictionary.
        Return value of post_kernel function is not required.
        '''
        pass
    
    def set(self, *args, **kwargs):
        '''
        set the attributes of the object.
        '''

        for k, v in kwargs.item():
            setattr(self, k, v)
        return self

    def run(self, *args, **kwargs):
        '''
        Call the kernel function of current object.  `args` will be passed
        to kernel function.  `kwargs` will be used to update the attributes of
        current object.  The return value of method run is the object itself.
        This allows a series of functions/methods to be executed in pipe.
        '''
        self.set(**kwargs)
        self.kernel(*args)
        return self

    __call__ = set

    def add_keys(self, **kwargs):
        '''Add or update attributes of the object and register these attributes in ._keys'''
        if kwargs:
            self.__dict__.update(**kwargs)
            self._keys = self._keys.union(kwargs.keys())
        return self

class Molecule(StreamObject):
    '''
    the basic class to hold 
    1) molecular structure
    2) gloabl options

    Attributes:
       verbose: init

       charge : int
       spin   : int or None
       atom : list
              format of atoms:
       unit : str: au [default]
       output : str: 
       basis : str ['AM1', 'PM6']

    '''

    def __init__(self, **kwargs):

        self.output = None

        self.verbose = 0
        self.charge = 0
        self.spin = 0 # 2j = n_alpha - n_beta

        self.atoms = [] 
        self.coordinates = None

        self.basis = None
        self.theory= 'AM1' # [DFT(NWchem), DFTB]

        return None


