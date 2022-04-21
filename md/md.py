"""
basic module for mds
Molecular dynamics (namd) module


"""
import torch

class md(torch.nn.Module):
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
