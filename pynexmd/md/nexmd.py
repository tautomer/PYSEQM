
#import torch
from .md import MD

class NEXMD(MD):
    def __init__(self, 
            mol,      # molecule 
            mdtype,   # TSH, MCE, AIMC, path-integral, etc.
            esdriver, # electronic structure driver
            *args,
            ):
        """
        esdriver: chose the driver for electronic structure:
                semi-empirical, DFT, CCSD (pySCF), QED-CCSD, HIPNN, etc.
        mdtype: "SH" (default), "AIMC", "Ehrenfest"
        args: dicitonary of other arguments: 
              1) classical time step
              2) thermostat
              3) number of classical steps
              4) number of excited states
              5) mdtypes 
              6) classical time step
              7) classic/quantum time ration
              8) decoherence type
              9) crossing check
              10) verbose
              ...
        """

        super().__init__()
        self.mdtype = mdtype


    def propagation(self):

        return None

    def check_crossing(self):
        return None

    def hopping(self):

        return None

    def output(self):

        return None
