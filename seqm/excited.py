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


def Vxi(mol, xi):

    return None

def tdhf_operator(mol, xi, method='rpa'):
    '''
    TDHF equations
    or RPA/CIS Liouville kernel

    corresponds to Lxi_testing subroutine in NEXMD fortran code
    
    xi Li(xi) == > L_ij
    Li = [F[rho0],x] + [V(xi),rho0]
             part1         part2

    part1 in MO (epsilon_a - epsilon_i) \delta_{iajb}
    Vxi in A.O.
    one-body in M.O.
    
    (A,   B) (x)  
    (-B, -A) (y)

    '''
    nocc = mol.nocc
    nvir = mol.nvir
    norb = mol.norb
    ncis = nocc * nvir

    #1) build Vxi in A.O.
    xi_ao = mo2ao(mol, xi)
    V_ao = Vxi(mol, xi_ao)
    
    # if solvent is required, add solvent (TODO)
    
    #2) tansform Vxi from AO TO MO
    # V^\dag A_AO V  --> A_MO

    V_mo = ao2mo(mol, V_ao)


    '''
    add one-body term 
    in MO representation: (epsilon_i - epsilon_a)\delta_{ij}\delta_{ab}
    '''
    i = 0
    for p in range(nocc):
        for h in range(nocc+1, norb):
            de  = mol.ehf[h] - mol.ehf[p]
            V_mol[i] += de * xi[i]
            V_mol[i+Ncis] = - V[i+Ncis] - de * xi[i]
            i += 1

    return V_mol

class TDHF():
    
    def kernel(mol, x0 = None):

        # set precondition
        precond = None

        if x0 is None:
            # set initial condition
            x0 = None (TODO)



        matvec = lambda vec: [tdhf_operator(mol, x) for x in vec]

        converged, evals, evecs = davidson(x0, precond, matvec,
                tol = self.conv_tol,
                nroots = nstates,
                maxcycle = self.maxcycle)
        
        return evals, evecs

def matvec_model(N,x):
    A = random matrix
    for i range(N):
       A[i,i] = 1.0
       A[i,i+1] = -1.0
       A[i+1,i] = -1.0

    return np.einsum('ij,j->i', A, x)


def Davidson(xi, precond, matvec, tol, nroots, maxcycle):

    
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
   
    """

    # ref: JCP xx

