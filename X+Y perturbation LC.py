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

def Pot_X3(N):
    V = 0
    #for k in range(1, N+1):
    V = sigmay_outer(1,N) * sigmax_k(1, N) + sigmay_outer(2,N) * sigmax_k(2, N) + sigmay_outer(3,N) * sigmax_k(3, N) 
    return V


def H_gad (V,N,L):
    H = H12_Ising(V,N) + L * (Pot_X(N)+Pot_Y(N))
    return H

def H_gad2 (V,N,L):
    H = H12_Ising(V,N) + L * (Pot_X(N))
    return H

def H_gad4 (V,N,L):
    H = H12_Ising(V,N) + L * (2*Pot_X(N))
    return H

def H_gad3 (V,N,L):
    H = H12_Ising(V,N) + L * (Pot_X(N)+Pot_X3(N))
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
    return (sigmax_comp(1,N) * sigmax_comp(2,N) * sigmax_comp(3,N)) \
        + (sigmay_comp(1,N) * sigmay_comp(2,N) * sigmay_comp(3,N))

def H_comp2 (N):
    return (sigmax_comp(1,N) * sigmax_comp(2,N) * sigmax_comp(3,N)) 

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

def H_right2 (N, L, k):
    return -(k * (-L)**k)/np.math.factorial(k-1) * tensor(H_comp2 (N), P_plus(N))


'#------------------------------------------------------------------------------'
# Creating effective Hamiltonian

def H_eff (d, L):
    H = 0.0
    es = H_gad(1/2, N_anc, L).eigenstates()
    for i in range (0,d):
        H += es[0][i] * (es[1][i] * es[1][i].dag())
    return H

def H_eff2 (d, L):
    H = 0.0
    es = H_gad2(1/2, N_anc, L).eigenstates()
    for i in range (0,d):
        H += es[0][i] * (es[1][i] * es[1][i].dag())
    return H

def H_eff3 (d, L):
    H = 0.0
    es = H_gad3(1/2, N_anc, L).eigenstates()
    for i in range (0,d):
        H += es[0][i] * (es[1][i] * es[1][i].dag())
    return H

def H_eff4 (d, L):
    H = 0.0
    es = H_gad4(1/2, N_anc, L).eigenstates()
    for i in range (0,d):
        H += es[0][i] * (es[1][i] * es[1][i].dag())
    return H

#defining energy-shift function and its projection

def Pi (d,L):
    P = 0.0
    es = H_gad(1/2, N_anc, L).eigenstates()
    for i in range (0,d):
        P += (es[1][i] * es[1][i].dag())
    return Qobj(np.array(P))

def f(L, k):
    a = 21
    b = 0.9
    c = 0.5
    d = 0
    return a*(L**k) + b*(L**(k-1)) + c*(L**(k-2)) + d

'#------------------------------------------------------------------------------'
#Find largest and smallest error

def DOM (A): #diagonal and ordered matrix
    Eigv = np.linalg.eig(A)
    sorted_eigv = np.sort(Eigv[0]) #smallest eigenvalue = first entry
    D = np.diag(sorted_eigv)
    return D

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
# Plot error of norm (compare paper p. 5) 
    
L1 = np.linspace(0.01, 0.2, 21)

def groundstate_diff (L):
    f0_vals = []
    f_vals = []
    for l in L:
        diff = Qobj((H_right(3, l, 3))).groundstate()[0] - Qobj(np.array(H_eff(2**3,l)) + np.array(f(l,3)*Pi(2**(3), l))).groundstate()[0]
        f_vals.append(np.abs(diff))
        diff2 = Qobj((H_right(3, l, 3))).groundstate()[0] - Qobj(np.array(H_eff(2**3,l))).groundstate()[0]
        f0_vals.append(np.abs(diff2))
    plt.plot(L, f0_vals, label = "f = 0")
    plt.plot(L, f_vals, label = "f = 21\u03BB$^3$ + 0.9\u03BB$^2$ + 0.5\u03BB")
    ax = plt.subplot(111)
    ax.spines['right'].set_color((.8,.8,.8))
    ax.spines['top'].set_color((.8,.8,.8))
    plt.title("X+Y Perturbation: Difference in groundstate energy", fontsize=16)
    plt.xlabel('Perturbation strength \u03BB', fontsize=14)
    plt.ylabel("$|E_0(H_{id}) - E_0(H_{eff})|$", fontsize=14)
    ax.grid(ls = "-.", c = "lightgrey")
    ax.legend()
    plt.savefig('X+Y groundstate-diff.png', dpi=300)

#groundstate_diff (L1)

def groundstate_separate (N, L):
    eff_val = []
    eff_val2 = []
    id_val = []
    for l in L:
        id_val.append(Qobj(np.array(H_right(N, l, 3))).groundstate()[0])
        eff_val.append(Qobj(np.array(H_eff(2**(N), l))).groundstate()[0])
        eff_val2.append(Qobj((np.array(H_eff(2**(N), l)) + np.array(f(l,3)*Pi(2**(3), l)))).groundstate()[0])
    plt.plot(L, id_val, color = "darkred", label = "H_id")
    plt.plot(L, eff_val, 'o', color = "grey", label = "H_eff, f = 0")
    plt.plot(L, eff_val2, 'o', color = "darkblue", label = "H_eff, f = 21\u03BB$^3$ + 0.9\u03BB$^2$ + 0.5\u03BB")
    ax = plt.subplot(111)
    ax.spines['right'].set_color((.8,.8,.8))
    ax.spines['top'].set_color((.8,.8,.8))
    plt.title("X+Y Perturbation: Groundstate energy comparison", fontsize=16)
    plt.xlabel('Perturbation strength \u03BB', fontsize=14)
    plt.ylabel("Groundstate energy $E_0$", fontsize=14)
    ax.grid(ls = "-.", c = "lightgrey")
    plt.legend()
    plt.savefig('X+Y groundstate-plot.png', dpi=300)
    
            
#groundstate_separate(3, L1)
 
def norm_diff (L):
    norm_vals = []
    norm_vals0 = []
    norm_vals2 = []
    norm_vals3 = []
    norm_vals4 = []
    for l in L:
        #diff = Qobj(np.array((H_right(3, l, 3)))) - (H_eff(2**3,l) + f(l,3)*Pi(2**(3), l))
        diff0 = Qobj(np.array((H_right(3, l, 3)))) - Qobj(np.array(H_eff(2**3,l)))
        diff2 = Qobj(np.array((H_right2(3, l, 3)))) - Qobj(np.array(H_eff2(2**3,l)))
        diff3 = Qobj(np.array((H_right(3, l, 3)))) - Qobj(np.array(H_eff3(2**3,l)))
        diff4 = Qobj(np.array((H_right(3, l, 3)))) - Qobj(np.array(H_eff4(2**3,l)))
        #norm_vals.append(op_norm(diff))
        norm_vals0.append(op_norm(diff0))
        norm_vals2.append(op_norm(diff2))
        norm_vals3.append(op_norm(diff3))
        norm_vals4.append(op_norm(diff4))
    #plt.plot(L, norm_vals, 'o', label = "f = 3.8 \u03BB$^3$ + 0.1 \u03BB$^2$")
    plt.plot(L, norm_vals0, label = "XXX+YYY ($H^{gad}_{X+Y}$)")
    plt.plot(L, norm_vals3, label = "XXX+YYY ($H^{gad}$)")
    plt.plot(L, norm_vals4, label = "XXX+XXX ($H^{gad}$)")
    plt.plot(L, norm_vals2, label = "XXX ($H^{gad}$)")
    ax = plt.subplot(111)
    ax.spines['right'].set_color((.8,.8,.8))
    ax.spines['top'].set_color((.8,.8,.8))
    plt.xlabel('Perturbation strength \u03BB', fontsize=14)
    plt.ylabel("$||H_{id} - H_{eff}||_2$", fontsize=14)
    ax.grid(ls = "-.", c = "lightgrey")
    plt.title("$H^{gad}$ vs $H^{gad}_{X+Y}$ : 2-norm error, f = 0", fontsize=16)
    ax.legend()     
    plt.savefig('X+Y 2-norm-diff.png', dpi=300)
    
norm_diff(L1)