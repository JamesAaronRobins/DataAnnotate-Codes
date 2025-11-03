'''Explanation of Code

The code starts with importing 2 modules. NumPy is used to define the bead coordinates as vectors, which simplifies the distance calculations. The second is JSON for the formatting of the final answer. 

Following this, the variables inputted in the prompt are added including the bead locations, ionic concentration, temperature and charge lengths. Below, more variables which contain various scientific constants are defined which include Avogadro constant, Boltzmann constant and unit conversions. 

The calculation starts with the temperature dependant dielectric constant which was given in the prompt before calculating the Bjerrum length and moving to the reduced charge (Q) given in the prompt as well. Below this the Debye length is calculated followed by the electrostatic perfactor (Std_Coef) calculation where units are also converted to kcal/mol for later in the calculation. 

After the prefactor the 3 different electrostatic coefficients are calculated for 1: the positive-negative interaction, 2: the positive-positive interaction and 3: the negative-negative interaction. The variable rcDist represents the inverse Deby Length (Kd) used in the Debye-Huckle equation. 

There are 6 pairwise distances calculated for each of the beads before the final electrostatic energy is calculated for each pair of beads and and summed up as ELEFinal. 

The last section contains the JSON formatting of the final answer with a new list "DistanceList" to hold each of the distances for the pairwise beads. 
'''
import numpy as np
import json 

#Variables
tempk = 298.15
negative_length_per_charge = 4.3817805
positive_length_per_charge = 2.8064008 
ionic_strength = 0.01
bead1 = np.array([34.5462, 67.4615, 66.5628])
bead2 = np.array([26.2351, 60.4741, 75.5781])
bead3 = np.array([44.6731, 72.6857, 71.7857])
bead4 = np.array([37.9672, 69.1341, 84.5721])
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

# Calculate reduced charge taking into account monovalent salt 
Zp = -(negative_length_per_charge / lb)#Negative Charges
posZP = 2*(positive_length_per_charge / lb)#Positive Charges

#Debye Length 
lambdaD = 1.0e10 * np.sqrt( (1.0e-3 * EPS0 * diele * BOLTZ_J) / (2.0 * N_AVO * ELE**2)  ) * np.sqrt(tempk / ionic_strength)
#Ele coefficient 
StdELE_coef = JOUL2KCAL_MOL * 1.0e10 * ELE**2 / (4.0e0 * np.pi * EPS0 * diele)

# Coefficients for each ELE intreaction
ele_coefmix = StdELE_coef*(posZP * Zp)
ele_coefpositive = StdELE_coef*(posZP * posZP)
ele_coefnegative = StdELE_coef*(Zp * Zp)

rcdist = 1.0 / lambdaD
e_ele = 0.0

#Pairwise distances for all beads 
dist12 = np.linalg.norm(bead1 - bead2)
dist13 = np.linalg.norm(bead1 - bead3)
dist14 = np.linalg.norm(bead1 - bead4)

dist23 = np.linalg.norm(bead2 - bead3)
dist24 = np.linalg.norm(bead2 - bead4)

dist34 = np.linalg.norm(bead3 - bead4)

#####################################################################################
ELEFinal = 0 # Defign the total ELE energy

#Ele Energy for beads 1 and 2 - Positive and Negative 
e_ele = (ele_coefmix/dist12*np.exp(-dist12*rcdist))
ELEFinal += e_ele
#print('1 ',e_ele)
#print('1 ',ELEFinal)

#Ele Energy for beads 1 and 3 - Positive and Negative 
e_ele = (ele_coefmix/dist13*np.exp(-dist13*rcdist))
ELEFinal += e_ele
#print('2 ',e_ele)
#print('2 ',ELEFinal)

#Ele Energy for beads 1 and 4 - Positive and Positive
e_ele = (ele_coefpositive/dist14*np.exp(-dist14*rcdist))
ELEFinal += e_ele
#print('3 ',e_ele)
#print('3 ',ELEFinal)

#Ele Energy for beads 2 and 3 - Negative and Negative
e_ele = (ele_coefnegative/dist23*np.exp(-dist23*rcdist))
ELEFinal += e_ele
#print('4 ',e_ele)
#print('4 ',ELEFinal)

#Ele Energy for beads 2 and 4 - Negative and Positive
e_ele = (ele_coefmix/dist24*np.exp(-dist24*rcdist))
ELEFinal += e_ele
#print(e_ele)
#print(ELEFinal)

#Ele Energy for beads 3 and 4 - Negative and Positive
e_ele = (ele_coefmix/dist34*np.exp(-dist34*rcdist))
ELEFinal += e_ele
#print(e_ele)
#print(ELEFinal)

DistanceList = [round(dist12,3), round(dist13,3), round(dist14,3), round(dist23,3), round(dist24,3), round(dist34,3)]
results = {
    "Dielectric constant": round(diele, 3),
    "Bjerrum length (nm)": round(lb/10, 3),
    "Charge on negative bead (e)": round(Zp, 3),
    "Charge on positive bead (e)": round(posZP, 3),
    "Debye length (nm)": round(lambdaD/10, 3),
    "Distance List (Angstrom) (1-2, 1-3, 1-4, 2-3, 2-4, 3-4)": DistanceList,
    "ELE Coefficient (Positive-Negative)": round(ele_coefmix,3),
    "ELE Coefficient (Positive-Positive)": round(ele_coefpositive,3),
    "ELE Coefficient (Negative-Negative)": round(ele_coefnegative,3),
    "Electrostatic energy (kcal/mol)": round(ELEFinal, 3)
}

# Print JSON-formatted output
print(json.dumps(results, indent=4))

'''Prompt

A system contains 4 electrostatically charged particles at positions: 
bead1 = 34.5462, 67.4615, 66.5628
bead2 = 26.2351, 60.4741, 75.5781
bead3 = 44.6731, 72.6857, 71.7857
bead4 = 37.9672, 69.1341, 84.5721
in a 0.01 M system. Bead 1 and bead 4 are positively charged with double the charge while beads 2 and 3 are negatively charged and all beads are in an implicit solvent at 298.15 K. What is the total electrostatic energy of the system in Kcal/mol? 
For reduced charges (Q) use the equation:
Q = \frac{b}{lB(T)}
Where b = 4.3817805 Angstrom for negative charges and b =  2.8064008 Angstrom for positive charges and lb(T) is the temperature dependant Bjerrum length. The temperature dependant dielectric constant can be calculated as: 
\epsilon(T) = 87.74-0.4008T_c+9.398*10^{-4}T_c^2-1.410*10^{-6}T_c^3
Where $T_c$ is the temperature in degrees Celsius. 
Include a step by step guide of all calculations to show the working. The final answer should be formatted as a JSON shown below with all numbers to 3 decimal places. The units for each measurement are shown next to their names:
Dielectric constant
Bjerrum length (nm)
Charge on negative bead (e)
Charge on positive bead (e)
Debye length (nm)
Distance List (Angstrom) (1-2, 1-3, 1-4, 2-3, 2-4, 3-4)
ELE Coefficient (Positive-Negative)
ELE Coefficient (Positive-Positive)
ELE Coefficient (Negative-Negative)
Electrostatic energy (kcal/mol)

'''
