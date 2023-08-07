#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 27 12:19:04 2023

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
             for i in range (2,N+k):
                 sig = tensor(sig,identity(2))
             sig = tensor(sig,sigmay())
             for i in range (k+1,N+1) :
                 sig = tensor(sig,identity(2))
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
# Creating effective Hamiltonian

def H_eff (d, L):
    H = 0.0
    es = H_gad(1/2, N_anc, L).eigenstates()
    for i in range (0,d):
        H += es[0][i] * (es[1][i] * es[1][i].dag())
    return Qobj(np.array(H))


#defining energy-shift function and its projection

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
#calculate the full Rydberg Hamiltonian

#define interatomic distances in µm
a_1 = 3.5
a_2 = 5
a_3 = 8.22
a_4 = 12.16

def H_Ryd_Ising (V_in_in, V_in_out1, V_in_out2, V_out_out, N=3):
    H_Ryd1 = 0.0  
    #inner qubit ising Iteraction:
    for k in range(1,N+1):  
        for j in range (k+1,N+1):  
           H_Ryd1 += V_in_in * (sigmaz_k(k,N) * sigmaz_k(j,N))

    #inner-(closest outer) Ising interaction
    for k in range(1,N+1):  
           H_Ryd1 += V_in_out1 * (sigmaz_k(k,N) * sigmaz_outer(k,N)) #i-ter mit (i+3)-ter Qubit
           
    #inner-(farest outer) Ising interaction
    for k in range(1, N+1):
        for j in range(1, N+1):
            if j == k:
                continue
            H_Ryd1 += V_in_out2 * (sigmaz_k(k,N) * sigmaz_outer(j,N))
    
    #outer-outer Ising interaction
    for k in range(N+1,N+1):  
        for j in range (k+1,N+1):  #4,5 - 4,6 - 5,6
            H_Ryd1 += V_out_out * (sigmaz_outer(k,N) * sigmaz_outer(j,N))        
           
    return H_Ryd1

def H_Ryd_trans (J_in_in, J_in_out1, J_in_out2, J_out_out, N=3):
    H_Ryd2 = 0.0  
    #inner qubit ising Iteraction:
    for k in range(1,N+1):  
        for j in range (k+1,N+1):  
           H_Ryd2 += J_in_in * (sigmax_k(k,N) * sigmax_k(j,N) + sigmay_k(k,N) * sigmay_k(j,N))

    #inner-(closest outer) Ising interaction
    for k in range(1,N+1):  
           H_Ryd2 += J_in_out1 * (sigmax_k(k,N) * sigmax_outer(j,N) + sigmay_k(k,N) * sigmay_outer(j,N)) #i-ter mit (i+3)-ter Qubit
           
    #inner-(farest outer) Ising interaction
    for k in range(1, N+1):
        for j in range(1, N+1):
            if j == k:
                continue
            H_Ryd2 += J_in_out2 * (sigmax_k(k,N) * sigmax_outer(j,N) + sigmay_k(k,N) * sigmay_outer(j,N))
    
    #outer-outer Ising interaction
    for k in range(N+1,N+1):  
        for j in range (k+1,N+1):  #4,5 - 4,6 - 5,6
            H_Ryd2 += J_out_out * (sigmax_outer(k,N) * sigmax_outer(j,N) + sigmay_outer(k,N) * sigmay_outer(j,N))      
           
    return H_Ryd2

def H_Ryd (V, J, N=3): #enter V,J as list of 4 values
    H_Ryd = H_Ryd_Ising(V[0], V[1], V[2], V[3]) + H_Ryd_trans(J[0], J[1], J[2], J[3])
    return H_Ryd

'#------------------------------------------------------------------------------'
#plot energy difference (error) between the two Hamiltonians

#let start with the example case from the powerpoint
V_list = [-2.26, 0, 0.007, 0.001]  #order: in-in| in-out1 | in-out2 | out-out
J_list = [0, 0.032, 0.007, 0.001]  #order: in-in| in-out1 | in-out2 | out-out
# define values for reduced Ryd-PTG (ONLY WANTED INTERACTIONS)
V0_list = [-2.26, 0, 0, 0]  #order: in-in| in-out1 | in-out2 | out-out
J0_list = [0, 0.032, 0, 0]  #order: in-in| in-out1 | in-out2 | out-out

    
L1 = np.linspace(0, 0.3, 31)

def Rydberg_vs_RydPTG(L):
    norm_vals = []
    norm_vals_PTG = []
    for l in L:
        # Update the J_list with the value of l
        J_list[1] = l
        J0_list[1] = l
        diff = Qobj(np.array(H_Ryd(V_list, J_list))) - Qobj(np.array(H_Ryd(V0_list, J0_list)))
        norm_vals.append(op_norm(diff))  
        x = -0.5/V_list[0]  #scaling factor for Ryd to PTG mapping
        diff_PTG = Qobj(np.array((H_right(3, l*x, 3)))) - (H_eff(2**3,l*x) + f(l*x,3)*Pi(2**(3), l*x))
        norm_vals_PTG.append(op_norm(diff_PTG)) 
    print(op_norm(diff))
    plt.plot(L, norm_vals, 'o', label ="Error between full Ryd Ham. and PTG-Ryd Ham")
    plt.plot(L, norm_vals_PTG, 'o', label = "Error between PTG-Ryd and 3-local Ham")
    ax = plt.subplot(111)
    ax.spines['right'].set_color((.8, .8, .8))
    ax.spines['top'].set_color((.8, .8, .8))
    plt.title("Full Rydberg vs PTG-Rydberg vs 3-local Hamiltonian")
    plt.xlabel("Inner-outer interaction strength")
    plt.ylabel("2-norm error")
    ax.grid(ls="-.", c="lightgrey")
    ax.legend()

Rydberg_vs_RydPTG(L1)

#map the values to the PTG values (V = -1/2):
# lambda = 0.007
l = 0.007
diff = Qobj(np.array((H_right(3, l, 3)))) - (H_eff(2**3,l) + f(l,3)*Pi(2**(3), l))
print(op_norm(diff))

def groundstate_diff (L):
    f_vals = []
    for l in L:
        diff = Qobj(np.array(H_eff(3,l)) + np.array(f(l,3)*Pi(2**(3), l))).groundstate()[0] - Qobj(np.array(H_Ryd(V_list, J_list))).groundstate()[0]
        f_vals.append(np.abs(diff))
    plt.plot(L, f_vals, label = "f = 4.1\u03BB$^3$ + 0.2\u03BB$^2$ + 0.1\u03BB$^2$")
    ax = plt.subplot(111)
    ax.spines['right'].set_color((.8,.8,.8))
    ax.spines['top'].set_color((.8,.8,.8))
    plt.title("Absolut energy difference between Ryd and PTG Hamiltonian")
    plt.xlabel('Perturbation strength \u03BB')
    plt.ylabel("Absolute error")
    ax.grid(ls = "-.", c = "lightgrey")
    ax.legend()

#groundstate_diff (L1)

def groundstate_separate (L):
    PTG_vals = []
    y2 = Qobj(np.array(H_Ryd(V_list, J_list))).groundstate()[0]
    for l in L:
        y1 = Qobj(np.array(H_eff(3,l)) + np.array(f(l,3)*Pi(2**(3), l))).groundstate()[0] 
        PTG_vals.append(np.abs(y1))
    plt.plot(L, PTG_vals, label = "PTG")
    plt.axhline(y=y2, color='red', linestyle='--', label='Ryd') 
    ax = plt.subplot(111)
    ax.spines['right'].set_color((.8,.8,.8))
    ax.spines['top'].set_color((.8,.8,.8))
    plt.title("Groudstate energy difference between Ryd and PTG Hamiltonian")
    plt.xlabel('Perturbation strength \u03BB')
    plt.ylabel("Groundstate energy")
    ax.grid(ls = "-.", c = "lightgrey")
    ax.legend()

#groundstate_separate(L1)

'#------------------------------------------------------------------------------'
#calculate and plot relative error
