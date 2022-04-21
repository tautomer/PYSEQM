import torch
from torch.autograd import grad
from .fock import fock
from .hcore import hcore
from .energy import elec_energy
from .SP2 import SP2
from .pack import *
from .diag import sym_eig_trunc, sym_eig_trunc1
import warnings
import time
#from .check import check
#scf_backward==0: ignore the gradient on density matrix
#scf_backward==1: use recursive formu
#scf_backward==2: go backward scf loop directly

debug=False
#debug=True

MAX_ITER = 1000
SCF_BACKWARD_MAX_ITER = 10
MAX_ITER_TO_STOP_IF_SCF_BACKWARD_DIVERGE = 5

RAISE_ERROR_IF_SCF_BACKWARD_FAILS = False
#if true, raise error rather than
#truncate gradient(set the gradient on non-converged molecules as 0.0)

RAISE_ERROR_IF_SCF_FORWARD_FAILS = False
#if true, raise error rather than ignore those non-convered molecules

def kernel():
    # kernel of TDHF equations

    return None


