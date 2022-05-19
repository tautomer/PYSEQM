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

def kernel(mol, xi):

    # kernel of TDHF equations
    # 
    # xi Li(xi) == > L_ij
    # Li = [F[rho0],x] + [V(xi),rho0]
    #
    # Vxi in A.O.
    # one-body in M.O.
    # 
    # (A,   B) (x)  
    # (-B, -A) (y)

    #1) build Vxi in A.O.

    #2) tansform Vxi from AO TO MO
    # V^\dag A_AO V  --> A_MO
    # 
    #) add one-body term

    return None


def Davidson(xi,matvec):

    
    """
    for m in range(k,mmax,k):
       if m <= k:
           for j in range(0,k):
               V[:,j] = t[:,j]/np.linalg.norm(t[:,j])
           theta_old = 1 
       elif m > k:
           theta_old = theta[:eig]
       V[:,:m],R = np.linalg.qr(V[:,:m])
       T = np.dot(V[:,:m].T, matvec(V[:,:m]))

       THETA,S = np.linalg.eig(T)
       idx = THETA.argsort()
       theta = THETA[idx]
       s = S[:,idx]
       for j in range(0,k):
           w = np.dot((A - theta[j]*I),np.dot(V[:,:m],s[:,j])) 
           q = w/(theta[j]-A[j,j])
           V[:,(m+j)] = q
       norm = np.linalg.norm(theta[:eig] - theta_old)
       if norm < tol:
           break
   
   #
       return None
    """



