#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  6 12:31:40 2023

@author: tillappel
"""
from arc import *
from IPython.display import display, HTML
import numpy as np
import scipy.constants as sc
import matplotlib.pyplot as plt

def find_largest_c3(n,n_2, l0, j0):
    largest_c3_d0 = 0
    largest_c3_d1 = 0
    largest_i_d0 = 0
    largest_i_d1 = 0
    largest_j_d0 = 0
    largest_j_d1 = 0
    largest_transition_d0 = ""
    largest_transition_d1 = ""

    atom = Rubidium()

    # Iterate over combinations of i and j
    for i in range(1, 4):
        for j in range(1, 4):
            # Calculate the dipole matrix element for pi/pi transition with d=0
            dsDME_pi_d0 = atom.getDipoleMatrixElement(n, l0, j0, j0, n+i, np.abs(l0-1), np.abs(j0-1), np.abs(j0), 0)
            dpDME_pi_d0 = atom.getDipoleMatrixElement(n_2, l0, j0, j0, n_2-j, l0+1, j0, j0, 0)
            c3_pi_d0 = (
                1
                / (4.0 * np.pi * sc.epsilon_0)
                * dsDME_pi_d0
                * dpDME_pi_d0
                * C_e**2
                * (sc.physical_constants["Bohr radius"][0]) ** 2
            )

            # Calculate the dipole matrix element for sigma+/sigma- transition with d=0
            dsDME_sigma_d0 = atom.getDipoleMatrixElement(n, l0, j0, j0, n+i, np.abs(l0-1), np.abs(j0-1), np.abs(j0),  -1)
            dpDME_sigma_d0 = atom.getDipoleMatrixElement(n_2, l0, j0, j0, n_2-j, l0+1, j0, j0, 1)

            c3_sigma_d0 = (
                1
                / (4.0 * np.pi * sc.epsilon_0)
                * dsDME_sigma_d0
                * dpDME_sigma_d0
                * C_e**2
                * (sc.physical_constants["Bohr radius"][0]) ** 2
            )

            # Compare the calculated c3 coefficients with d=0 and update the largest values
            if abs(c3_pi_d0) > abs(largest_c3_d0):
                largest_c3_d0 = c3_pi_d0
                largest_i_d0 = i
                largest_j_d0 = j
                largest_transition_d0 = "pi/pi"

            if abs(c3_sigma_d0) > abs(largest_c3_d0):
                largest_c3_d0 = c3_sigma_d0
                largest_i_d0 = i
                largest_j_d0 = j
                largest_transition_d0 = "sigma+/sigma-"

            # Calculate the dipole matrix element for pi/pi transition with d=1
            dsDME_pi_d1 = atom.getDipoleMatrixElement(n, l0, j0, j0, n+i, np.abs(l0-1), np.abs(j0-1), np.abs(j0-1), 0)
            dpDME_pi_d1 = atom.getDipoleMatrixElement(n_2, l0, j0, j0, n_2-j, l0+1, j0+1, j0+1, 0)

            c3_pi_d1 = (
                1
                / (4.0 * np.pi * sc.epsilon_0)
                * dsDME_pi_d1
                * dpDME_pi_d1
                * C_e**2
                * (sc.physical_constants["Bohr radius"][0]) ** 2
            )

            # Calculate the dipole matrix element for sigma+/sigma- transition with d=1
            dsDME_sigma_d1 = atom.getDipoleMatrixElement(n, l0, j0, j0, n+i, np.abs(l0-1), np.abs(-1+j0), np.abs(-1+j0),  -1)
            dpDME_sigma_d1 = atom.getDipoleMatrixElement(n_2, l0, j0, j0, n_2-j, l0+1, 1+j0, 1+j0, 1)

            c3_sigma_d1 = (
                1
                / (4.0 * np.pi * sc.epsilon_0)
                * dsDME_sigma_d1
                * dpDME_sigma_d1
                * C_e**2
                * (sc.physical_constants["Bohr radius"][0]) ** 2
            )

            # Compare the calculated c3 coefficients with d=1 and update the largest values
            if abs(c3_pi_d1) > abs(largest_c3_d1):
                largest_c3_d1 = c3_pi_d1
                largest_i_d1 = i
                largest_j_d1 = j
                largest_transition_d1 = "pi/pi"

            if abs(c3_sigma_d1) > abs(largest_c3_d1):
                largest_c3_d1 = c3_sigma_d1
                largest_i_d1 = i
                largest_j_d1 = j
                largest_transition_d1 = "sigma+/sigma-"

    return (
        largest_i_d0, largest_j_d0, largest_transition_d0, abs(largest_c3_d0) / C_h * 1.0e9,
        largest_i_d1, largest_j_d1, largest_transition_d1, abs(largest_c3_d1) / C_h * 1.0e9
    )


# Specify the value of n, l0, and j0
n = 59
n_2 = 59
l = 0
j = 0.5

# Find the largest C3 coefficients for d=0 and d=1, and their corresponding i, j, and transition
largest_i_d0, largest_j_d0, largest_transition_d0, largest_c3_d0, largest_i_d1, largest_j_d1, largest_transition_d1, largest_c3_d1 = find_largest_c3(n, n_2, l, j)

# Print the results
print("For d=0:")
print("Largest C3 of Rb(%dP -> %dS/%dD) = %.3f GHz (µm)^3 (i = %d, j = %d, Transition = %s)" % (n, n-largest_i_d0, n+largest_j_d0, largest_c3_d0, largest_i_d0, largest_j_d0, largest_transition_d0))

print("For d=1:")
print("Largest C3 of Rb(%dP -> %dS/%dD) = %.3f GHz (µm)^3 (i = %d, j = %d, Transition = %s)" % (n, n-largest_i_d1, n+largest_j_d1, largest_c3_d1, largest_i_d1, largest_j_d1, largest_transition_d1))

'--------------------------------------------------'
#resonant interaction of groundstate to excited state with opposite parity

atom = Rubidium(cpp_numerov=False)
dme = atom.getDipoleMatrixElement(63, 1, 1/2, 1/2, 40, 0, 1/2, 1/2, +1)

c3_2 = (
    1
    / (4.0 * np.pi * sc.epsilon_0)
    * dme
    * dme
    * C_e**2
    * (sc.physical_constants["Bohr radius"][0]) ** 2
)
print("C_3 of Rb(63 S -> 61P) = %.3f GHz (mu m)^3 " % (abs(c3_2) / C_h * 1.0e9))

'================================================='

# Evaluation of the Cs 60S_1/2 C6 coefficient using perturbation theory (Theta=0,phi=0)
l0 = 0
j0 = 0.5
mj0 = 0.5
# Target State
theta = 0
# Polar Angle [0-pi]
phi = 0
# Azimuthal Angle [0-2pi]
dn = 5
# Range of n to consider (n0-dn:n0+dn)
deltaMax = 25e9  # Max pair-state energy difference [Hz]

# Set target-state and extract value
calculation = PairStateInteractions(
    Rubidium(), n, l0, j0, n, l0, j0, mj0, mj0
)
C6 = calculation.getC6perturbatively(theta, phi, dn, deltaMax)
print("C6 [%s] = %.2f GHz (mum)^6" % (printStateString(n, l0, j0), C6))


'--------------------------------------------------'

# Define a range of values for n
n_values = range(30, 80)
a_1 = 1 #µm

# Lists to store the C3 and C6 coefficients for d=0 and d=1
c3_values_d0 = []
c3_values_d1 = []
c6_values = []

# Iterate over the values of n
for n in n_values:
    # Find the largest C3 coefficients for d=0 and d=1, and their corresponding i, j, and transition
    largest_i_d0, largest_j_d0, largest_transition_d0, largest_c3_d0, largest_i_d1, largest_j_d1, largest_transition_d1, largest_c3_d1 = find_largest_c3(n, n_2, l0, j0)

    # Append the largest C3 coefficients to the respective c3_values lists
    c3_values_d0.append(largest_c3_d0  / a_1**3)
    c3_values_d1.append(largest_c3_d1  / a_1**3)

    # Calculate the C6 coefficient
    calculation = PairStateInteractions(
        Rubidium(), n, l0, j0, n, l0, j0, mj0, mj0
    )
    C6 = calculation.getC6perturbatively(theta, phi, dn, deltaMax)

    # Append the C6 coefficient to the c6_values list
    c6_values.append(np.abs(C6) / a_1**6)

#Plotting the C3 and C6 coefficientsplt.plot(n_values, c3_values_d1, label="Largest C3 Coefficient")
#plt.plot(n_values, c3_values_d1, label="C3 Coefficient (d=1)")
#plt.plot(n_values, c6_values, label="C6 Coefficient")
'-------------------'
plt.semilogy(n_values, c3_values_d0, label="Largest C3 Coefficient") #CURRENTLY: d=1
plt.semilogy(n_values, c6_values, label="C6 Coefficient")
'-------------------'
plt.xlabel("n")
plt.ylabel("C3, C6 [GHz]")
plt.legend(fontsize = "large", loc="upper left")
plt.title("C3 & C6 coefficients of Rb |n,S>")
plt.savefig('log plot S c3,c6.png', dpi=300)
plt.show()


