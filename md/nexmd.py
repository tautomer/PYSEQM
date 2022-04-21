"""
basic module for mds
Molecular dynamics (namd) module


"""
import torch
from .md import md

class nexmd(md):
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

        super().__init__()


    def propagation(self):

        return None

    def check_crossing(self):
        return None

    def hopping(self):

        return None

    def output(self):

        return None
