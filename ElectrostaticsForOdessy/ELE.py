
import numpy as np
import json 

#Variables
tempk = 298.15
length_per_charge = 4.3817805
polychargelength = 2.8064008 
ionic_strength = 0.01
bead1 = np.array([34.5462, 67.46153, 66.5628])
bead2 = np.array([26.2351, 60.47415, 75.5781])
#Fixed Constants
MM_A=87.740e0               #1
MM_B=-0.4008e0              #2
MM_C=9.398e-4               #3
MM_D=-1.410e-6              #4 - These 1-4 are the constants from "Dielectric constant of water from 0 to 100 C, Malmberg & Maryott, 1956"
EPS0 = 8.8541878128e-12   # Electric constant [F/m] (Vaccume permiativity)
ELE  = 1.602176634e-19    # Elementary charge [C]
BOLTZ_J = 1.380649e-23    # Boltzmann constant [J/K]
N_AVO   = 6.02214076e23   # Avogadro constant [/mol]
KCAL2JOUL = 4184.0              #(kcal -> J)  [J/kcal]
JOUL2KCAL = 1.0/KCAL2JOUL   #< (J -> kcal)  [kcal/J]
JOUL2KCAL_MOL  = JOUL2KCAL * N_AVO  #< (J -> kcal/mol)
BOLTZ_KCAL_MOL = BOLTZ_J * JOUL2KCAL_MOL   #< Boltzmann constant [kcal/mol/K]

#Start Calculation
# Temperature Dependant dielectric consstant  
Tc = tempk - 273.15
diele =  MM_A + MM_B*Tc + MM_C*Tc*Tc + MM_D*Tc*Tc*Tc

#Bjerrum length
lb = ELE * ELE / (4.0e0 * np.pi * EPS0 * diele * BOLTZ_J * tempk) * 1.0e10

# Calculate reduced charge taking into account monovalent salt condensation
#xi = lb / length_per_charge
#theta = 1.0 - 1.0/xi
Zp = -(length_per_charge / lb)#RNA is negative 
polyZP = 2*(polychargelength / lb)#PolymerisPositive

#Debye Length 
lambdaD = 1.0e10 * np.sqrt( (1.0e-3 * EPS0 * diele * BOLTZ_J) / (2.0 * N_AVO * ELE**2)  ) * np.sqrt(tempk / ionic_strength)
#Ele coefficient 
StdELE_coef = JOUL2KCAL_MOL * 1.0e10 * ELE**2 / (4.0e0 * np.pi * EPS0 * diele)
# Coefficient between bead 1 and 2
ele_coef = StdELE_coef*(polyZP * Zp)

rcdist = 1.0 / lambdaD
e_ele = 0.0
ELEInterchain = 0.0

dist = np.linalg.norm(bead1 - bead2)

e_ele = e_ele + (ele_coef/dist*np.exp(-dist*rcdist))

Eele = e_ele

#print('Dielectric constant (H2O): ', diele)
#print("LB = ", lb)
#print("Debye Length", lambdaD)
#print("Negative Charge", Zp)
#print("Positive Charge", polyZP)
#print("Distance", dist)
#print("FINAL ENERGY = ", round(Eele,3))

results = {
    "Dielectric constant": round(diele, 3),
    "Bjerrum length (nm)": round(lb/10, 3),
    "Charge on negative bead (e)": round(Zp, 3),
    "Charge on positive bead (e)": round(polyZP, 3),
    "Debye length (nm)": round(lambdaD/10, 3),
    "Distance between beads (nm)": round(dist/10, 3),
    "Electrostatic energy (kcal/mol)": round(Eele, 3)
}

# Print JSON-formatted output
print(json.dumps(results, indent=4))

'''print('##########################################')
print('CHATGPT CODE BELOW')
print('##########################################')

from math import pi, exp, sqrt
# Constants
T_K = 298.15  # Kelvin
T_C = T_K - 273.15  # Celsius
salt_conc_M = 0.01  # 10 mM
e = 1.602176634e-19  # Elementary charge (C)
epsilon_0 = 8.854187817e-12  # Vacuum permittivity (F/m)
Na = 6.02214076e23  # Avogadro's number
k_B = 1.380649e-23  # Boltzmann constant (J/K)

# Length per charge in Ångström, convert to meters
b_neg_A = 4.3817805
b_pos_A = 2.8064008
b_neg = b_neg_A * 1e-10
b_pos = b_pos_A * 1e-10

# Bead positions in Ångström
r1 = np.array([34.5462, 67.46153, 66.5628])
r2 = np.array([26.2351, 60.47415, 75.5781])
r_ij_m = np.linalg.norm(r1 - r2) * 1e-10  # Distance in meters

# Dielectric constant ε(T)
epsilon_T = (
    87.74
    - 0.4008 * T_C
    + 9.398e-4 * T_C**2
    - 1.410e-6 * T_C**3
)

# Bjerrum length in meters
lB = e**2 / (4 * pi * epsilon_0 * epsilon_T * k_B * T_K)

# Charges (dimensionless, in units of e)
Q_neg = - (b_neg / lB)'''
'''Important Change - CHAPTGPT GAVE POSITIVE VALUE FOR NEG BEAD - DIDNT RECOGNISE THE CHARGE INVERSION STEP'''

'''Q_pos = 2 * (b_pos / lB)  # Positive bead is 2x stronger

# Debye length in meters
I = salt_conc_M *1000 # Ionic strength for 1:1 electrolyte '''
'''Important Change - CHAPTGPT GAVE INCORRECT I VALUE - IT was 1000 OFF DUE TO INCORRECT UNITS'''
'''lambda_D = sqrt(epsilon_T * epsilon_0 * k_B * T_K / (2 * Na * e**2 * I))
K_D = 1 / lambda_D

# Electrostatic energy in Joules
q1 = Q_neg * e
q2 = Q_pos * e
U_J = (q1 * q2) / (4 * pi * epsilon_0 * epsilon_T) * exp(-K_D * r_ij_m) / r_ij_m
#print("DEBUG - ",U_J*Na)
# Convert to kcal/mol
U_kcalmol = U_J * Na / 4184

# Print results
print(f"Dielectric constant ε(T): {epsilon_T:.3f}")
print(f"Bjerrum length: {lB * 1e9:.3f} nm")
print(f"Charge on negative bead: {Q_neg:.3f} e")
print(f"Charge on positive bead: {Q_pos:.3f} e")
print(f"Debye length: {lambda_D * 1e9:.3f} nm")
print(f"Distance between beads: {r_ij_m * 1e9:.3f} nm")
print(f"Electrostatic energy: {U_kcalmol:.3f} kcal/mol")

'''