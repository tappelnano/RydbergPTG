#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 30 13:48:20 2023

@author: tillappel

"""

from qutip import *
import numpy as np
import scipy as sc
import matplotlib.pyplot as plt
from scipy.optimize import minimize_scalar


'#-----------------------------------------------------------------------------'
#General definitions and functions


def sigmaz_k(k,N): #Pauli-Z Interaction of inner (ancilla) qubits
    if k<=N:
         if k==1:
             sig = sigmaz()
             for i in range (2,2*N+1):  #counting includes outer qubits
                 sig = tensor(sig,identity(2))
         elif k==N :
             sig = identity(2)
             for i in range (2,N):
                 sig = tensor(sig,identity(2))
             sig = tensor(sig,sigmaz())
             for j in range(2,N+2):
                 sig = tensor(sig,identity(2))
         else :
             sig = identity(2)
             for i in range (2,k):
                 sig = tensor(sig,identity(2))
             sig = tensor(sig,sigmaz())
             for i in range (k+1,2*N+1) :
                 sig = tensor(sig,identity(2))
         return sig
    else :  
         print("k can't be greater than N, wrong input dude")
         
        
def sigmax_k(k,N): 
    if k<=N:
         if k==1:
             sig = sigmax()
             for i in range (2,2*N+1):
                 sig = tensor(sig,identity(2))
                 
         elif k==N :
             sig = identity(2)
             for i in range (2,N):
                 sig = tensor(sig,identity(2))
             sig = tensor(sig,sigmax())
             for j in range (2,N+2):
                 sig = tensor(sig,identity(2))
         else :
             sig = identity(2)
             for i in range (2,k):
                 sig = tensor(sig,identity(2))
             sig = tensor(sig,sigmax())
             for i in range (k+1,2*N+1) :
                 sig = tensor(sig,identity(2))
         return sig
    else :  
         print("k can't be greater than N, wrong input dude")
   
def sigmay_k(k,N): 
    if k<=N:
         if k==1:
             sig = sigmay()
             for i in range (2,2*N+1):
                 sig = tensor(sig,identity(2))
                 
         elif k==N :
             sig = identity(2)
             for i in range (2,N):
                 sig = tensor(sig,identity(2))
             sig = tensor(sig,sigmay())
             for j in range (2,N+2):
                 sig = tensor(sig,identity(2))
         else :
             sig = identity(2)
             for i in range (2,k):
                 sig = tensor(sig,identity(2))
             sig = tensor(sig,sigmay())
             for i in range (k+1,2*N+1) :
                 sig = tensor(sig,identity(2))
         return sig
    else :  
         print("k can't be greater than N, wrong input dude")         
   
def sigmaz_outer(k,N): 
    if k<=N:
         if k==1:
             sig = sigmaz()
             for i in range (2,N+2):
                 sig = tensor(identity(2),sig)
             for j in range (2,N+1):
                 sig = tensor(sig, identity(2))  #überprüfen!!
         elif k==N :
             sig = identity(2)
             for i in range (2,2*N):
                 sig = tensor(sig,identity(2))
             sig = tensor(sig,sigmaz())
         else :
             sig = identity(2)
             for i in range (2,N+k):
                 sig = tensor(sig,identity(2))
             sig = tensor(sig,sigmaz())
             for i in range (k+1,N+1) :
                 sig = tensor(sig,identity(2))
         return sig
    else :  
         print("k can't be greater than N, wrong input dude")         

def sigmax_outer(k,N): 
    if k<=N:
         if k==1:
             sig = sigmax()
             for i in range (2,N+2):
                 sig = tensor(identity(2),sig)
             for j in range (2,N+1):
                 sig = tensor(sig, identity(2))  #überprüfen!!
         elif k==N :
             sig = identity(2)
             for i in range (2,2*N):
                 sig = tensor(sig,identity(2))
             sig = tensor(sig,sigmax())
         else :
             sig = identity(2)
             for i in range (2,N+k):
                 sig = tensor(sig,identity(2))
             sig = tensor(sig,sigmax())
             for i in range (k+1,N+1) :
                 sig = tensor(sig,identity(2))
         return sig
    else :  
         print("k can't be greater than N, wrong input dude")  

def sigmay_outer(k,N): 
    if k<=N:
         if k==1:
             sig = sigmay()
             for i in range (2,N+2):
                 sig = tensor(identity(2),sig)
             for j in range (2,N+1):
                 sig = tensor(sig, identity(2))  #überprüfen!!
         elif k==N :
             sig = identity(2)
             for i in range (2,2*N):
                 sig = tensor(sig,identity(2))
             sig = tensor(sig,sigmay())
         else :
             sig = identity(2)
             for i in range (2,N+k):
                 sig = tensor(sig,identity(2))
             sig = tensor(sig,sigmay())
             for i in range (k+1,N+1) :
                 sig = tensor(sig,identity(2))
         return sig
    else :  
         print("k can't be greater than N, wrong input dude")  
         
def IdeN(N):  
     sig = identity(2)  
     for i in range (2,N+1) :  
         sig = tensor(sig,identity(2))  
     return sig  

def H12_Ising(V,N):      #modified Ising Hamiltonian (set J = 1/2) 
     H12 = 0.0  
     for k in range(1,N+1):  
         for j in range (k+1,N+1):  
            H12 += IdeN(2*N) - (sigmaz_k(k,N) * sigmaz_k(j,N))
     return V*H12

def Pot_X(N):
    V = 0
    #for k in range(1, N+1):
    V = sigmax_outer(1,N) * sigmax_k(1, N) + sigmax_outer(2,N) * sigmax_k(2, N) + sigmax_outer(3,N) * sigmax_k(3, N) 
    return V

def Pot_Y(N):
    V = 0
    #for k in range(1, N+1):
    V = sigmay_outer(1,N) * sigmay_k(1, N) + sigmay_outer(2,N) * sigmay_k(2, N) + sigmay_outer(3,N) * sigmay_k(3, N) 
    return V

def H_gad (V,N,L):
    H = H12_Ising(V,N) + L/2 * (Pot_X(N)+Pot_Y(N))
    return H

#Set of Variables / System-Parameters
N_anc = 3   #number of inner (or outer) qubits
lam = 0.1    #perturbation strenght (parameter lambda)

'#------------------------------------------------------------------------------'
# H_comp and corresponding operators

def sigmaz_comp(k,N): #Pauli-Z Interaction of inner (ancilla) qubits
    if k<=N:
         if k==1:
             sig = sigmaz()
             for i in range (2,N+1):  #counting includes outer qubits
                 sig = tensor(sig,identity(2))
         elif k==N :
             sig = identity(2)
             for i in range (2,N):
                 sig = tensor(sig,identity(2))
             sig = tensor(sig,sigmaz())
         else :
             sig = identity(2)
             for i in range (2,k):
                 sig = tensor(sig,identity(2))
             sig = tensor(sig,sigmaz())
             for i in range (k+1,N+1) :
                 sig = tensor(sig,identity(2))
         return sig
    else :  
         print("k can't be greater than N, wrong input habub")
             
def sigmax_comp(k,N): 
    if k<=N:
         if k==1:
             sig = sigmax()
             for i in range (2,N+1):
                 sig = tensor(sig,identity(2))
         elif k==N :
             sig = identity(2)
             for i in range (2,N):
                 sig = tensor(sig,identity(2))
             sig = tensor(sig,sigmax())
         else :
             sig = identity(2)
             for i in range (2,k):
                 sig = tensor(sig,identity(2))
             sig = tensor(sig,sigmax())
             for i in range (k+1,N+1) :
                 sig = tensor(sig,identity(2))
         return sig
    else :  
         print("k can't be greater than N, wrong input habub x")
   
def sigmay_comp(k,N): 
    if k<=N:
         if k==1:
             sig = sigmay()
             for i in range (2,N+1):
                 sig = tensor(sig,identity(2))
                 
         elif k==N :
             sig = identity(2)
             for i in range (2,N):
                 sig = tensor(sig,identity(2))
             sig = tensor(sig,sigmay())
         else :
             sig = identity(2)
             for i in range (2,k):
                 sig = tensor(sig,identity(2))
             sig = tensor(sig,sigmay())
             for i in range (k+1,N+1) :
                 sig = tensor(sig,identity(2))
         return sig
    else :  
         print("k can't be greater than N, wrong input habub y")  

def H_comp1 (N):
    return 0.5*(sigmax_comp(1,N) * sigmax_comp(2,N) * sigmax_comp(3,N)) \
        + 0.5*(sigmay_comp(1,N) * sigmay_comp(2,N) * sigmay_comp(3,N))


'#------------------------------------------------------------------------------'
# P+ Projector and right side of the theorem

def P_plus (N):
    ket_1 = basis(2)
    for i in range (1, N):
        ket_1 = tensor(ket_1, basis(2))
    ket_0 = basis(2,1)
    for i in range (1, N):
        ket_0 = tensor(ket_0, basis(2,1))
    plus = 1/np.sqrt(2) * (ket_1 + ket_0)
    P = plus * plus.dag()
    return P

def H_right (N, L, k):
    return -(k * (-L)**k)/np.math.factorial(k-1) * tensor(H_comp1 (N), P_plus(N))

'#------------------------------------------------------------------------------'
# Creating effective Hamiltonian, Pi (Support-projector), and the energy shift function f(lambda)

def H_eff (d, L):
    H = 0.0
    es = H_gad(1/2, N_anc, L).eigenstates()
    for i in range (0,d):
        H += es[0][i] * (es[1][i] * es[1][i].dag())
    return Qobj(np.array(H))

def Pi (d,L):
    P = 0.0
    es = H_gad(1/2, N_anc, L).eigenstates()
    for i in range (0,d):
        P += (es[1][i] * es[1][i].dag())
    return Qobj(np.array(P))

def f(L, k):
    a = 4.1
    b = 0.2
    c = 0.1
    d = 0
    return a*(L**k) + b*(L**(k-1)) + c*(L**(k-2)) + d

def op_norm (A):
    max_val = 0
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
#Find largest and smallest error

def DOM (A): #diagonal and ordered matrix
    Eigv = np.linalg.eig(A)
    sorted_eigv = np.sort(Eigv[0]) #smallest eigenvalue = first entry
    D = np.diag(sorted_eigv)
    return D

#example for lambda = .1
LEFT = H_eff(2**(3),0.1) 
RIGHT = H_right(3, 0.1, 3)
DIFF = DOM(LEFT)-DOM(RIGHT)

LEFT_F2 = H_eff(2**(3),0.1) + f(0.1,3)*Pi(2**(3), 0.1)
RIGHT = H_right(3, 0.1, 3)
DIFF_F2 = DOM(LEFT_F2)-DOM(RIGHT)

gs_diff = DIFF[0][0] #difference in groundstate energy
max_diff = abs(DIFF).max()

print(f"The groundstate energy difference is: {gs_diff}")
print(f"The largest energy difference is: {max_diff}")

def error_source(matrices, labels = None):
    n = max(len(matrix) for matrix in matrices)  # Get the maximum size of the matrices

    # Plotting
    for i, matrix in enumerate(matrices):
       y = np.diag(matrix) # Extract the diagonal elements
       x = range(1, len(y) + 1)  # Generate x-axis values
       if labels is not None:
           label = labels[i]
       else:
           label = f'Matrix {i+1}'
       plt.plot(x, y, marker='o', label= label)
    plt.xlabel('Excitation level')
    plt.ylabel('Energy difference')
    plt.title('H_eff - H_id, d = 8')
    plt.legend()
    plt.show()

#error_source([DIFF, DIFF_F2], ['f = 0', 'f = 4.1 \u03BB$^3$ + 0.2 \u03BB$^2$ + 0.1 \u03BB'])

'#------------------------------------------------------------------------------'
# Plot error of norm (compare paper p. 5)
    
L1 = np.linspace(0.01, 0.2, 21)

def rel_error_plot (L):
    norm_vals = []
    norm_vals0 = []
    for l in L:
        diff = Qobj(np.array((H_right(3, l, 3)))) - (H_eff(2**3,l) + f(l,3)*Pi(2**(3), l))
        diff0 = Qobj(np.array((H_right(3, l, 3)))) - (H_eff(2**3,l))
        norm_vals.append(op_norm(diff))
        norm_vals0.append(op_norm(diff0))
    plt.plot(L, norm_vals, label = "f = 4.1\u03BB$^3$ + 0.2\u03BB$^2$ + 0.1\u03BB")
    plt.plot(L, norm_vals0, label = "f = 0")
    ax = plt.subplot(111)
    ax.spines['right'].set_color((.8,.8,.8))
    ax.spines['top'].set_color((.8,.8,.8))
    plt.xlabel('Perturbation strength \u03BB', fontsize=14)
    plt.title("X+Y Perturbation: 2-norm error", fontsize=16)
    plt.ylabel("$||H_{id} - H_{eff}||_2$", fontsize=14)
    ax.grid(ls = "-.", c = "lightgrey")
    ax.legend()     
    plt.savefig('X+Y 2-norm-diff.png', dpi=300)
    
rel_error_plot (L1)

def rel_error_bounds (L):
    upper_bound = []
    lower_bound = []
    for l in L:
        diff = Qobj(np.array((H_right(3, l, 3)))) - (H_eff(2**3,l) + f(l,3)*Pi(2**(3), l))
        upper_bound.append(op_norm(diff))
        lower_bound.append(np.linalg.norm(diff, ord = -2))
    plt.plot(L, (upper_bound), label = "$||H_{id} - H_{eff}||_2$", color = "darkred")
    plt.plot(L, (lower_bound), label = "$||H_{id} - H_{eff}||_{-2}$", color = "darkgreen")
    ax = plt.subplot(111)
    ax.spines['right'].set_color((.8,.8,.8))
    ax.spines['top'].set_color((.8,.8,.8))
    plt.title("X+Y Perturbation: 2-norm error bounds", fontsize=16)
    plt.xlabel('Perturbation strength \u03BB', fontsize=14)
    plt.ylabel("($\pm$)2-norm error", fontsize=14)
    ax.grid(ls = "-.", c = "lightgrey")
    ax.fill_between(L, upper_bound, lower_bound, color='green', alpha=.2)
    ax.legend()
    plt.savefig('X+Y error bounds.png', dpi=300)
    

#rel_error_bounds(L1)
