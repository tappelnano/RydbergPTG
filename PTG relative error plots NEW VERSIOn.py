#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct  8 12:24:47 2023

@author: tillappel
"""

from qutip import *
import numpy as np
import scipy as sc
import matplotlib.pyplot as plt
from scipy.optimize import minimize_scalar

'#-----------------------------------------------------------------------------'
# Definition of Pauli-Operators

# (1) Pauli-Operators for inner/ancilla qubits (H_anc) (2^2Nx2^2N dim):
    
def sigmax_anc (k,N): #N = number of ancilla qubits
    if k <= N:
        if k==1:
             sig = sigmax()
             for i in range (1,2*N): #runs til 2N = number of ancilla+comp. qubits
                 sig = tensor(sig,identity(2)) #sigmax x I x ... x I
        elif k==N:
             sig = identity(2)
             for i in range (1,N-1):
                 sig = tensor(sig,identity(2))
             sig = tensor(sig,sigmax())
             for j in range (1,N+1):
                 sig = tensor(sig,identity(2))
        else :
             sig = identity(2)
             for i in range (2,k):
                 sig = tensor(sig,identity(2))
             sig = tensor(sig,sigmax())
             for i in range (k,2*N) :
                 sig = tensor(sig,identity(2))  
        return sig #dimension output matric = 2^(2N) x 2^(2N)
    else :  
         print("k can't be greater than N") #error checker

def sigmay_anc(k,N): #same construction as simgax_anc
    if k<=N:
         if k==1:
             sig = sigmay()
             for i in range (1,2*N):
                 sig = tensor(sig,identity(2))
                 
         elif k==N :
             sig = identity(2)
             for i in range (1,N-1):
                 sig = tensor(sig,identity(2))
             sig = tensor(sig,sigmay())
             for j in range (1,N+1):
                 sig = tensor(sig,identity(2))
         else :
             sig = identity(2)
             for i in range (2,k):
                 sig = tensor(sig,identity(2))
             sig = tensor(sig,sigmay())
             for i in range (k,2*N) :
                 sig = tensor(sig,identity(2))
         return sig
    else :  
         print("k can't be greater than N") 
         
def sigmaz_anc(k,N): 
    if k<=N:
         if k==1:
             sig = sigmaz()
             for i in range (1,2*N):
                 sig = tensor(sig,identity(2))
         elif k==N :
             sig = identity(2)
             for i in range (1,N-1):
                 sig = tensor(sig,identity(2))
             sig = tensor(sig,sigmaz())
             for j in range(1,N+1):
                 sig = tensor(sig,identity(2))
         else :
             sig = identity(2)
             for i in range (2,k):
                 sig = tensor(sig,identity(2))
             sig = tensor(sig,sigmaz())
             for i in range (k,2*N) :
                 sig = tensor(sig,identity(2))
         return sig
    else :  
         print("k can't be greater than N, wrong input dude")
         
#(2) Pauli-Operators for computational/outer qubits (V) (2^2Nx2^2N dim)

def sigmax_outer(k,N): 
    if k<=N:
         if k==1:
             sig = sigmax()
             for i in range (1,N+1):
                 sig = tensor(identity(2),sig)
             for j in range (1,N):
                 sig = tensor(sig, identity(2))  #überprüfen!!
         elif k==N :
             sig = identity(2)
             for i in range (1,2*N-1):
                 sig = tensor(sig,identity(2))
             sig = tensor(sig,sigmax())
         else :
             sig = identity(2)
             for i in range (1,N+k-1):
                 sig = tensor(sig,identity(2))
             sig = tensor(sig,sigmax())
             for i in range (k,N) :
                 sig = tensor(sig,identity(2))
         return sig
    else :  
         print("k can't be greater than N")  

def sigmay_outer(k,N): 
    if k<=N:
         if k==1:
             sig = sigmay()
             for i in range (1,N+1):
                 sig = tensor(identity(2),sig)
             for j in range (1,N):
                 sig = tensor(sig, identity(2))  #überprüfen!!
         elif k==N :
             sig = identity(2)
             for i in range (1,2*N-1):
                 sig = tensor(sig,identity(2))
             sig = tensor(sig,sigmay())
         else :
             sig = identity(2)
             for i in range (1,N+k-1):
                 sig = tensor(sig,identity(2))
             sig = tensor(sig,sigmay())
             for i in range (k,N) :
                 sig = tensor(sig,identity(2))
         return sig
    else :  
         print("k can't be greater than N") 
         
def sigmaz_outer(k,N): 
    if k<=N:
         if k==1:
             sig = sigmaz()
             for i in range (1,N+1):
                 sig = tensor(identity(2),sig)
             for j in range (1,N):
                 sig = tensor(sig, identity(2))  #überprüfen!!
         elif k==N :
             sig = identity(2)
             for i in range (1,2*N-1):
                 sig = tensor(sig,identity(2))
             sig = tensor(sig,sigmaz())
         else :
             sig = identity(2)
             for i in range (2,N+k):
                 sig = tensor(sig,identity(2))
             sig = tensor(sig,sigmaz())
             for i in range (k,N) :
                 sig = tensor(sig,identity(2))
         return sig
    else :  
         print("k can't be greater than N")      
         
#(3) Pauli-Operators for H_comp (2^Nx2^N dim.):
    
def sigmax_comp(k,N): 
    if k<=N:
         if k==1:
             sig = sigmax()
             for i in range (1,N):
                 sig = tensor(sig,identity(2))
         elif k==N :
             sig = identity(2)
             for i in range (1,N-1):
                 sig = tensor(sig,identity(2))
             sig = tensor(sig,sigmax())
         else :
             sig = identity(2)
             for i in range (1,k-1):
                 sig = tensor(sig,identity(2))
             sig = tensor(sig,sigmax())
             for i in range (k,N) :
                 sig = tensor(sig,identity(2))
         return sig
    else :  
         print("k can't be greater than N")
   
def sigmay_comp(k,N): 
    if k<=N:
         if k==1:
             sig = sigmay()
             for i in range (1,N):
                 sig = tensor(sig,identity(2))
                 
         elif k==N :
             sig = identity(2)
             for i in range (1,N-1):
                 sig = tensor(sig,identity(2))
             sig = tensor(sig,sigmay())
         else :
             sig = identity(2)
             for i in range (1,k-1):
                 sig = tensor(sig,identity(2))
             sig = tensor(sig,sigmay())
             for i in range (k,N) :
                 sig = tensor(sig,identity(2))
         return sig
    else :  
         print("k can't be greater than N")  
         
def sigmaz_comp(k,N): 
    if k<=N:
         if k==1:
             sig = sigmaz()
             for i in range (1,N):  
                 sig = tensor(sig,identity(2))
         elif k==N :
             sig = identity(2)
             for i in range (1,N-1):
                 sig = tensor(sig,identity(2))
             sig = tensor(sig,sigmaz())
         else :
             sig = identity(2)
             for i in range (1,k-1):
                 sig = tensor(sig,identity(2))
             sig = tensor(sig,sigmaz())
             for i in range (k,N) :
                 sig = tensor(sig,identity(2))
         return sig
    else :  
         print("k can't be greater than N")

'#-----------------------------------------------------------------------------'
#H_gad = H_anc + V operators:

def H_anc(N):      #modified Ising Hamiltonian (J = -1/2) for ancilla qubits
     H = 0.0  
     for k in range(1,N+1):  
         for j in range (k+1,N+1):  
            H += identity(2**(2*N)) - np.array((sigmaz_anc(k,N) * sigmaz_anc(j,N)))
     return (1/2)*H
 
def V_XYZ(N): 
    V = sigmax_outer(1,N) * sigmax_anc(1, N) + sigmax_outer(2,N) * sigmay_anc(2, N) + sigmax_outer(3,N) * sigmaz_anc(3, N) 
    return np.array(V)

def H_gad (N,L):
    H = H_anc(N) + L * V_XYZ(N)
    return H
    
'#-----------------------------------------------------------------------------'
#H_comp and RHS of the equation

def H_comp_XYZ (N):
    H = sigmax_comp(1,N) * sigmay_comp(2,N) * sigmaz_comp(3,N)
    return np.array(H)

def P_plus (N): #acts on *ancilla* register
    ket_1 = basis(2)
    for i in range (1, N):
        ket_1 = tensor(ket_1, basis(2))
    ket_0 = basis(2,1)
    for i in range (1, N):
        ket_0 = tensor(ket_0, basis(2,1))
    plus = 1/np.sqrt(2) * (ket_1 + ket_0)
    P = plus * plus.dag()
    return np.array(P)

def H_id (N, L, k):
    H = -(k * (-L)**k)/np.math.factorial(k-1) * tensor(Qobj(H_comp_XYZ (N)), Qobj(P_plus (N)))
    return np.array(H)

'#-----------------------------------------------------------------------------'
#Defining the effective Hamiltonian with the support-projector Pi and f(lambda)

def H_eff (N, d, L): #order of spectral decomposition
    H = 0.0
    es = H_gad(N, L).eigenstates()
    for i in range (0,d):
        H += es[0][i] * (es[1][i] * es[1][i].dag())
    return np.array(H)

def Pi (N, d,L): #projection onto supports (d-th order)
    P = 0.0
    es = H_gad(N, L).eigenstates()
    for i in range (0,d):
        P += (es[1][i] * es[1][i].dag())
    return np.array(P)

def f(L, k): #f(lambda)
    a = 0
    b = 0.995
    c = 0.995
    d = 0
    return a*(L**k) + b*(L**(k-1)) + c*(L**(k-2)) + d

def f2(L, k): #f(lambda)
    a = 0
    b = 1.05
    c = 1
    d = 0
    return a*(L**k) + b*(L**(k-1)) + c*(L**(k-2)) + d

def f3(L, k): #f(lambda)
    a = 0
    b = 1
    c = 1
    d = 0
    return a*(L**k) + b*(L**(k-1)) + c*(L**(k-2)) + d

def op_norm (M):
    max_val = 0
    A = Qobj(M)
    gs = A.groundstate()[1]
    hs = (-A).groundstate()[1]
    val1 = np.absolute(expect(A, hs))
    val2 = np.absolute(expect(A, gs))
    if val2 > val1:
        max_val = val2
    else:    
        max_val = val1
    return max_val

'#------------------------------------------------------------------------------'
# Plot error of norm (compare paper p. 5)
    
L1 = np.linspace(0.01, 0.15, 21)

def rel_error_plot_XYZ (L, d): #L = list with perturbation strengths (-> lambda)
    norm_vals = []
    norm_vals2 = []
    norm_vals3 = []
    for l in L:
        diff = H_id(3, l, 3) - (H_eff(3, d,l) + f(l,3)*Pi(3, d, l))  #chose d=2*3=8 for the q ancilla qubit case
        diff2 = H_id(3, l, 3) - (H_eff(3, d,l) + f2(l,3)*Pi(3, d, l))
        diff3 = H_id(3, l, 3) - (H_eff(3, d,l) + f3(l,3)*Pi(3, d, l))
        norm_vals.append(   op_norm(diff) / op_norm(H_id(3,l,3))  )
        norm_vals2.append(   op_norm(diff2) / op_norm(H_id(3,l,3))  )
        norm_vals3.append(   op_norm(diff3) / op_norm(H_id(3,l,3))  )
    plt.plot(L, norm_vals, label = "f  = 0.995 \u03BB$^2$ + 0.995 \u03BB")
    plt.plot(L, norm_vals2, label = "f = 1.05 \u03BB$^2$ + \u03BB")
    plt.plot(L, norm_vals3, label = "f = \u03BB$^2$ + \u03BB")
    ax = plt.subplot(111)
    ax.spines['right'].set_color((.8,.8,.8))
    ax.spines['top'].set_color((.8,.8,.8))
    plt.title("XYZ Perturbation: 2-norm error", fontsize=16)
    plt.xlabel('Perturbation strength \u03BB', fontsize=14)
    plt.ylabel("$\\frac{||H_{id} - H_{eff}||_2}{||H_{id}||}$", fontsize=14)
    ax.grid(ls = "-.", c = "lightgrey")
    ax.legend()
    
rel_error_plot_XYZ (L1, 2**3) #chose d=2*3=8 for the q ancilla qubit case