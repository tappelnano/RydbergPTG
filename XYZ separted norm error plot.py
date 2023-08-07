#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  9 18:23:45 2023

@author: tillappel
"""

from qutip import *
import numpy as np
import matplotlib.pyplot as plt

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

def Pot_ext(N):
    V = 0
    #for k in range(1, N+1):
    V = sigmax_outer(1,N) * sigmax_k(1, N) + sigmay_outer(2,N) * sigmax_k(2, N) + sigmay_outer(3,N) * sigmax_k(3, N)
    return V

def H_gad (V,N,L):
    H = H12_Ising(V,N) + L * Pot_ext(N)
    return H

#Definition of measurements
def z_Basis (state, N):    #Pauli-Z measurement
    psi = []
    for n in range (1, N+1):
        psi.append(expect(Qobj(np.array(sigmaz_k(n,N))), Qobj(np.array(state))))
    return psi

def x_Basis (state, N):   #Pauli-X measurement
    psi = []
    for n in range (1, N+1):
        psi.append(expect(Qobj(np.array(sigmax_k(n,N))), Qobj(np.array(state))))
    return psi

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
         print("k can't be greater than N, wrong input dude")
         
        
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
         print("k can't be greater than N, wrong input dude")
   
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
         print("k can't be greater than N, wrong input dude")  

def H_comp4 (N):
    return sigmax_comp(1,N) * sigmay_comp(2,N) * sigmaz_comp(3,N)

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
    return -(k * (-L)**k)/np.math.factorial(k-1) * tensor(H_comp4 (N), P_plus(N))

'#------------------------------------------------------------------------------'
# Creating effective Hamiltonian

def H_eff (d, L):
    H = 0.0
    es = H_gad(1/2, N_anc, L).eigenstates()
    for i in range (0,d):
        H += es[0][i] * (es[1][i] * es[1][i].dag())
    return H
# paper says d = 2**N, with N = number of qubits on which H_comp is operating
fock_dm = np.array(H_eff(2**(N_anc), 0.1))     #2**3 = 8 
dm_es = Qobj(np.around(fock_dm, 4)).eigenstates()

'#------------------------------------------------------------------------------'
# Plot error of norm (compare paper p. 5)
    
L1 = np.linspace(0, 0.15, 11)

def error_separate (N, L):
    for l in L:
        H_ef = H_eff(2**N, l)
        H_id = H_right(N, l, 3)
        y = H_id.norm(norm = 'tr') - H_ef.norm(norm = 'tr')  
        plt.plot(l, -y, 'o', color = "darkblue")
        ax = plt.subplot(111)
        ax.spines['right'].set_color((.8,.8,.8))
        ax.spines['top'].set_color((.8,.8,.8))
        plt.xlabel('Perturbation strength \u03BB')
        plt.ylabel("Error")
        ax.grid(ls = "-.", c = "lightgrey")

