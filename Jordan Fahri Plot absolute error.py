
"""
Created on Wed Mar  1 16:12:15 2023

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

def Pot_XYZ(N):
    V = 0
    #for k in range(1, N+1):
    V = sigmax_outer(1,N) * sigmax_k(1, N) + sigmay_outer(2,N) * sigmax_k(2, N) + sigmaz_outer(3,N) * sigmax_k(3, N)
    return V  #Qobj(np.array(V))

def Pot_XYZ_XYY(N):
    V = 0
    #for k in range(1, N+1):
    V = sigmax_outer(1,N) * sigmax_k(1, N) + sigmay_outer(2,N) * sigmax_k(2, N) + sigmaz_outer(3,N) * sigmax_k(3, N) \
        + sigmax_outer(1,N) * sigmax_k(1, N) + sigmay_outer(2,N) * sigmax_k(2, N) + sigmay_outer(3,N) * sigmax_k(3, N) 
    return V  #Qobj(np.array(V))



def Pot_XYZZ(N):
    V = 0
    #for k in range(1, N+1):
    V = sigmax_outer(1,N) * sigmax_k(1, N) + sigmay_outer(2,N) * sigmax_k(2, N) + sigmaz_outer(3,N) * sigmax_k(3, N) \
        + sigmaz_outer(4,N) * sigmax_k(4, N) 
    return V  #Qobj(np.array(V))

def H_gad_XYZ (V,N,L):
    H = H12_Ising(V,N) + L * Pot_XYZ(N)
    return Qobj(np.array(H))

def H_gad_XYZ_XYY (V,N,L):
    H = H12_Ising(V,N) + L * Pot_XYZ_XYY(N)
    return Qobj(np.array(H))

def H_gad_XYZZ (V,N,L):
    H = H12_Ising(V,N) + L * Pot_XYZZ(N)
    return Qobj(np.array(H))

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

def H_comp_XYZ (N):
    return sigmax_comp(1,N) * sigmay_comp(2,N) * sigmaz_comp(3,N)

def H_comp_XYZ_XYY (N):
    return sigmax_comp(1,N) * sigmay_comp(2,N) * sigmaz_comp(3,N) \
        + sigmax_comp(1,N) * sigmay_comp(2,N) * sigmay_comp(3,N)

def H_comp_XYZZ (N):
    return sigmax_comp(1,N) * sigmay_comp(2,N) * sigmaz_comp(3,N)* sigmaz_comp(4,N)

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
    
def H_right_XYZ (N, L, k):
    return -(k * (-L)**k)/np.math.factorial(k-1) * tensor(H_comp_XYZ (N), P_plus(N))

def H_right_XYZ_XYY (N, L, k):
    return -(k * (-L)**k)/np.math.factorial(k-1) * tensor(H_comp_XYZ_XYY (N), P_plus(N))

def H_right_XYZZ (N, L, k):
    return -(k * (-L)**k)/np.math.factorial(k-1) * tensor(H_comp_XYZZ (N), P_plus(N))

'#------------------------------------------------------------------------------'
# Creating effective Hamiltonian, Pi (Support-projector), and the energy shift function f(lambda)

def H_eff_XYZ (d, L):
    H = 0.0
    es = H_gad_XYZ(1/2, 3, L).eigenstates()
    for i in range (0,d):
        H += es[0][i] * (es[1][i] * es[1][i].dag())
    return Qobj(np.array(H))

def H_eff_XYZ_XYY (d, L):
    H = 0.0
    es = H_gad_XYZ_XYY(1/2, 3, L).eigenstates()
    for i in range (0,d):
        H += es[0][i] * (es[1][i] * es[1][i].dag())
    return Qobj(np.array(H))

def H_eff_XYZZ (d, L):
    H = 0.0
    es = H_gad_XYZZ(1/2, 4, L).eigenstates()
    for i in range (0,d):
        H += es[0][i] * (es[1][i] * es[1][i].dag())
    return Qobj(np.array(H))

def Pi_XYZ (d,L):
    P = 0.0
    es = H_gad_XYZ(1/2, N_anc, L).eigenstates()
    for i in range (0,d):
        P += (es[1][i] * es[1][i].dag())
    return Qobj(np.array(P))

def Pi_XYZ_XYY (d,L):
    P = 0.0
    es = H_gad_XYZ_XYY(1/2, N_anc, L).eigenstates()
    for i in range (0,d):
        P += (es[1][i] * es[1][i].dag())
    return Qobj(np.array(P))

def Pi_XYZZ (d,L):
    P = 0.0
    es = H_gad_XYZZ(1/2, 4, L).eigenstates()
    for i in range (0,d):
        P += (es[1][i] * es[1][i].dag())
    return Qobj(np.array(P))

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

def f(L, k):
    a = 4.1# 151.5
    b = 0.2
    c = 0.1
    d = 0
    return a*(L**k) + b*(L**(k-1)) + c*(L**(k-2)) + d


'#------------------------------------------------------------------------------'
# Plot error of norm (compare paper p. 5)
    
L1 = np.linspace(0.01, 0.15, 21)
A = np.linspace(0, -3, 101)

def error_plot_f0 (L):
    norm_valsXYZ = []
    norm_valsXYZ_XYY = []
    norm_valsXYZZ = []
    for l in L:
        diffXYZ = Qobj(np.array((H_right_XYZ(3, l, 3)))) - (H_eff_XYZ(2**3,l))
        diffXYZ_XYY = Qobj(np.array((H_right_XYZ_XYY(3, l, 3)))) - (H_eff_XYZ_XYY(2**3,l))
        diffXYZZ = Qobj(np.array((H_right_XYZZ(4, l, 4)))) - (H_eff_XYZZ(2**4,l))
        norm_valsXYZ.append(op_norm(diffXYZ)) 
        norm_valsXYZ_XYY.append(op_norm(diffXYZ_XYY)) 
        norm_valsXYZZ.append(op_norm(diffXYZZ)) 
    plt.plot(L, (norm_valsXYZ), label = "XYZ")
    plt.plot(L, (norm_valsXYZ_XYY), label = "XYZ + XYY")
    plt.plot(L, (norm_valsXYZZ), label = "XYZZ")
    ax = plt.subplot(111)
    ax.spines['right'].set_color((.8,.8,.8))
    ax.spines['top'].set_color((.8,.8,.8))
    plt.xlabel('Perturbation strength \u03BB', fontsize=14)
    plt.ylabel("$||H_{id} - H_{eff}||_2$", fontsize=14)
    plt.title("PTG Plots: 2-norm, f = 0", fontsize=16)
    ax.grid(ls = "-.", c = "lightgrey")
    ax.legend()
    plt.savefig('PTG plot f0.png', dpi=300)
    
#error_plot_f0(L1)    
    
def error_plot_f (L):
    norm_valsXYZ = []
    norm_valsXYZ_XYY = []
    norm_valsXYZZ = []
    for l in L:
        diffXYZ = Qobj(np.array((H_right_XYZ(3, l, 3)))) - (H_eff_XYZ(2**3,l) + f(l,3)*Pi_XYZ(2**(3), l))
        diffXYZ_XYY = Qobj(np.array((H_right_XYZ_XYY(3, l, 3)))) - (H_eff_XYZ_XYY(2**3,l) + f(l,3)*Pi_XYZ_XYY(2**(3), l))
        diffXYZZ = Qobj(np.array((H_right_XYZZ(4, l, 4)))) - (H_eff_XYZZ(2**4,l) + f(l,4)*Pi_XYZZ(2**(4), l))
        norm_valsXYZ.append(op_norm(diffXYZ)) 
        norm_valsXYZ_XYY.append(op_norm(diffXYZ_XYY)) 
        norm_valsXYZZ.append(op_norm(diffXYZZ)) 
    plt.plot(L, (norm_valsXYZ), label = "XYZ")
    plt.plot(L, (norm_valsXYZ_XYY), label = "XYZ + XYY")
    plt.plot(L, (norm_valsXYZZ), label = "XYZZ")
    ax = plt.subplot(111)
    ax.spines['right'].set_color((.8,.8,.8))
    ax.spines['top'].set_color((.8,.8,.8))
    plt.xlabel('Perturbation strength \u03BB', fontsize=14)
    plt.ylabel("$||H_{id} - H_{eff}||_2$", fontsize=14)
    plt.title("PTG Plots: 2-norm, f = 4.1\u03BB$^3$ + 0.2\u03BB$^2$ + 0.1\u03BB", fontsize=16)
    ax.grid(ls = "-.", c = "lightgrey")
    ax.legend()
    plt.savefig('PTG plot wf.png', dpi=300)
    
#error_plot_f(L1)      
    
def error_plot_rel (L):
    norm_valsXYZ = []
    norm_valsXYZ_XYY = []
    norm_valsXYZZ = []
    for l in L:
        diffXYZ = Qobj(np.array((H_right_XYZ(3, l, 3)) / l**3)) - (H_eff_XYZ(2**3,l))
        diffXYZ_XYY = Qobj(np.array((H_right_XYZ_XYY(3, l, 3))/ l**3 )) - (H_eff_XYZ_XYY(2**3,l))
        diffXYZZ = Qobj(np.array((H_right_XYZZ(4, l, 4))/ l**4 )) - (H_eff_XYZZ(2**4,l))
        norm_valsXYZ.append(op_norm(diffXYZ)) 
        norm_valsXYZ_XYY.append(op_norm(diffXYZ_XYY)) 
        norm_valsXYZZ.append(op_norm(diffXYZZ)) 
    plt.plot(L, (norm_valsXYZ), label = "XYZ")
    plt.plot(L, (norm_valsXYZ_XYY), label = "XYZ + XYY")
    plt.plot(L, (norm_valsXYZZ), label = "XYZZ")
    ax = plt.subplot(111)
    ax.spines['right'].set_color((.8,.8,.8))
    ax.spines['top'].set_color((.8,.8,.8))
    plt.xlabel('Perturbation strength \u03BB', fontsize=14)
    plt.ylabel("$||H_{id} - H_{eff}||_2$", fontsize=14)
    plt.title("PTG Plots: 2-norm, f = 0", fontsize=16)
    ax.grid(ls = "-.", c = "lightgrey")
    ax.legend()
    plt.savefig('PTG plot rel.png', dpi=300)
    

error_plot_rel(L1)

