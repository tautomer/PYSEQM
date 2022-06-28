
import os, sys
import torch
#import bml
import numpy as np

#TODO
def dot_product(A,B):

    '''
    solver: torch, bml, numpy
    datataype: dense, ellpack, csr, etc
    dp: single, double, 
    '''

    C = None

    return C

def inner_product(A,B):

    C = None

    return C


def davidson(x0, precond, matvec, tol=1.e-7, nroots=3, max_cycle=200, max_space=100, verbose=1):

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
   
    converged = False
    evals = None
    evecs = None

    return converged, evals, evecs


 
if __name__ == '__main__':

    N = 1000
    # generate a symmetric matrix A
    # and use davidson to calculate the eigenvalues and eigenvectors of A
    nroots = 5
    A = np.zeros((N,N))
    for i in range(N):
        A[i,i] = 0.1
        if i < N-1:
            A[i,i+1] = 1.0
            A[i+1,i] = 1.0


    x0 = None # TODO
    precond = None

    matvec = lambda vec: [np.einsum('ij,j->i', A, x) for x in vec]

    converged, e, v = davidson(x0, precond, matvec)
