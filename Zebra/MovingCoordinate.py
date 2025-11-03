''' 
Code Explanation 

This code begins with importing 2 modules, numpy and copy. Numpy is used for calculations and copy is used to make copies deep copies of lists. The next line of code, np.random.seed(27), is used to ensure the random number seed is always the same, for repeatability when running the code multiple times. 

The coordinates of bead1 and bead2 are then defined as well as copies of these used at the end of the code for calculating the difference between the inital and final positions. Numerous variables are then defined for steps, dt, friction (fric), mass 1, mass 2 and temperature in Kelvin. The constants are then defined with the force constant (k) and distance value (r0) before the physical constants and unit conversions such as Avogadro and Boltzmann constants. More variables are labelled below, though only to set them before the loop starts. 

The functions used for calculation are below. A function RND_Bump is defined to generate random numbers, but for this specific problem, it is commented out. Instead, the specific bump values provided in the prompt are used to ensure a reproducible result. The second function sets the MD_Coef values by implementing the outline from the prompt of the MD coefficients function. This is used for the Langevin integration steps. This function is hard coded for 2 particles for simplicity of calling the function, though can be edited to be called as many times as needed with simple modification. The MD_coefficients are made of 4 parts: 
md_coef[0]: Scales the unitless random number (bump) into a physical stochastic acceleration. It incorporates mass, temperature, friction, and the timestep.
md_coef[1]: A dimensionless factor that damps the velocity from the previous step due to friction. If friction were zero, this value would be 1.
md_coef[2]: Converts the deterministic force (in kcal/mol/angstrom) into an acceleration component for the velocity update.
md_coef[3]: Used in the second half of the velocity-Verlet algorithm to update the positions (x_new = x_old + v_new * md_coef[3]).

The force calculation for the harmonic bond is next which is derived from the energy equation given in the prompt. This function calculates the stochastic acceleration (the "random kick") for each bead. It uses md_coef[0] to scale the provided unitless random bump values into physically meaningful acceleration terms. The final function labelled as velocity is the velocity-Verlet-style Langevin integration (Brünger-Brooks-Karplus (BBK) integrator) where the 3 sections (old velocity, deterministic force, and random force) which are friction damped previous velocities, the deterministic forces, and and a stochastic term that averages the random acceleration from the previous and current step. 

The loop to determine coordinates then starts. This begins with a for loop ranging to 'steps', a variable set at the top of the script. The random bumps are first defined. These were given in the prompt ensuring repeatability. The MD_coef values are then set for the langevin integration, followed by calculating acceleration (random force) and then the deterministic forces (bond). Finally, the velocity is called for the Langevin integration. 

After the functions are called, deep copies of the velocities and accelerations are saved. They are deep copies to create new copies of the lists, rather than regular copies to avoid editing the lists when the values are recalculated in future loops (not needed here as only 1 step is taken, but good practice). The final lines calculate the movement of the coordinates in the bead1diff and bead2diff sections before adding this difference to the coordinates for the final coordinate values. 

The last section outside the loop is the formatting for the printing of the answer. 

'''
import numpy as np
import copy

np.random.seed(27)

# Bead Positions
bead1 = np.array([743.2363, 45.6731, 153.5673])
bead2 = np.array([741.0813, 48.6721, 158.6731])
bead1start = np.array([743.2363, 45.6731, 153.5673])
bead2start = np.array([741.0813, 48.6721, 158.6731])
#Variables 
Steps = 1
dt = 0.2
fric = 0.0001
mass1 = 328.212
mass2 = 344.212
tempk = 298.15

#Constants for Force Equation
k = 15.0
r0 = 6.13

#Physics Constants
BOLTZ_J = 1.380649e-23    # Boltzmann constant [J/K]
N_AVO   = 6.02214076e23   # Avogadro constant [/mol]
KCAL2JOUL = 4184.0              #(kcal -> J)  [J/kcal]
JOUL2KCAL = 1.0/KCAL2JOUL   #< (J -> kcal)  [kcal/J]
JOUL2KCAL_MOL  = JOUL2KCAL * N_AVO  #< (J -> kcal/mol)
BOLTZ_KCAL_MOL = BOLTZ_J * JOUL2KCAL_MOL   #< Boltzmann constant [kcal/mol/K]
#print(BOLTZ_KCAL_MOL)

# Updates for each step
force = [0,0,0]
Forces1 = 0
Forces2 = 0
XYZ_Move = 0
prevvelo1 = [0,0,0]
prevvelo2 = [0,0,0]
prevacceleration1 = [0,0,0]
prevacceleration2 = [0,0,0]
# Functions for Calculation

#Random Bump - Sets the random numbers for the random bump - These are also given in the prompt and are not neeeded to be recalculated
# however for further use of the script (e.g. different systems) the random numbers can be calculated again
def RND_Bump():
    x = np.random.normal()
    y = np.random.normal()
    z = np.random.normal()
    rnd_bump = [x,y,z]
    return rnd_bump

# The molecular dynamics coefficients for the langevin dynamics integration - Both MDCoef 1 and 2 are calculated in the same coefficient
# This is hard coded to just be 2 particles but is easily modified for larger systems
def set_MD_Coef():
    #MD Coef for bead 1
    md_coef1 = [0,0,0,0]
    c1 = 0.5 * dt * fric / mass1
    c2 = np.sqrt(1.0 / (1.0 + c1))

    #md_coef(1) = sqrt(b) / 2m * sqrt(2 gamma kT h)
    md_coef1[0] = 0.5 * c2 / mass1 * np.sqrt(2.0 * fric * BOLTZ_KCAL_MOL * tempk * dt)
    #md_coef(2) = a
    md_coef1[1] = (1.0 - c1) / (1.0 + c1)
    #md_coef(3) = sqrt(b) h / m
    md_coef1[2] = c2 * dt / mass1
    #md_coef(4) = sqrt(b) h
    md_coef1[3] = c2 * dt
    
    #MD Coef for bead 2
    md_coef2 = [0,0,0,0]
    c1 = 0.5 * dt * fric / mass2
    c2 = np.sqrt(1.0 / (1.0 + c1))

    #md_coef(1) = sqrt(b) / 2m * sqrt(2 gamma kT h)
    md_coef2[0] = 0.5 * c2 / mass2 * np.sqrt(2.0 * fric * BOLTZ_KCAL_MOL * tempk * dt)
    #md_coef(2) = a
    md_coef2[1] = (1.0 - c1) / (1.0 + c1)
    #md_coef(3) = sqrt(b) h / m
    md_coef2[2] = c2 * dt / mass2
    #md_coef(4) = sqrt(b) h
    md_coef2[3] = c2 * dt
    
    return md_coef1, md_coef2

# Calculate the bond force from the harmoic potential given in the prompt
def Calc_Force(bead1, bead2):
    dist12 = np.linalg.norm(bead1 - bead2)
    vector = bead1 - bead2
    #print(vector)
    delta = dist12-r0
    force = k * (delta/dist12) * vector
    Forces1 = -force
    Forces2 = +force
    return Forces1, Forces2

# Randoom force / bump used in langevin integration 
def Accelerations(bump1, bump2, md_coef1, md_coef2):
    accels_pre1 = [0,0,0]
    accels_pre2 = [0,0,0]
    # bump is a dimensionless value which needs to be scaled to a physical force. This is done here
    accels_pre1[0] = md_coef1[0] * bump1[0]
    accels_pre1[1] = md_coef1[0] * bump1[1]
    accels_pre1[2] = md_coef1[0] * bump1[2]
    
    accels_pre2[0] = md_coef2[0] * bump2[0]
    accels_pre2[1] = md_coef2[0] * bump2[1]
    accels_pre2[2] = md_coef2[0] * bump2[2]
    
    return accels_pre1, accels_pre2

# Langevin integration of the velocity 
def Calc_Velocity(prevvelo, prevacceleration, md_coef, forces, acceleration):
    #velocity = ( old velocity) + (deterministic force) + ( random force)
    velocity = [0,0,0]
    velocity[0] =  md_coef[1] * prevvelo[0] + md_coef[2] * forces[0] + (prevacceleration[0] + acceleration[0])
    velocity[1] =  md_coef[1] * prevvelo[1] + md_coef[2] * forces[1] + (prevacceleration[1] + acceleration[1])
    velocity[2] =  md_coef[1] * prevvelo[2] + md_coef[2] * forces[2] + (prevacceleration[2] + acceleration[2])
    return velocity

# Calculate the force felt by each bead based on current positions - This loop can be extended for a number of steps (e.g. 1, 1000, 1000000)
for i in range(Steps):
    
    #bump1 = RND_Bump()
    #bump2 = RND_Bump()
    bump1 = [1.2856054152322811, -0.30355337730615073, 0.6190756618694475]
    bump2 = [0.39599854914377247, 0.22340564775539556, -0.05433941722861124]
    
    md_coef1, md_coef2 = set_MD_Coef()

    acceleration1, acceleration2 = Accelerations(bump1,bump2, md_coef1,md_coef2)
    
    forces1, forces2 = Calc_Force(bead1, bead2)

    velocity1 = Calc_Velocity(prevvelo1, prevacceleration1, md_coef1, forces1, acceleration1)
    velocity2 = Calc_Velocity(prevvelo2, prevacceleration2, md_coef2, forces2, acceleration2)
        
    prevvelo1 = copy.deepcopy(velocity1)
    prevvelo2 = copy.deepcopy(velocity2)
    prevacceleration1 = copy.deepcopy(acceleration1)
    prevacceleration2 = copy.deepcopy(acceleration2)
    
    bead1diff = [0,0,0]
    bead1diff[0] = md_coef1[3]*velocity1[0]
    bead1diff[1] = md_coef1[3]*velocity1[1]
    bead1diff[2] = md_coef1[3]*velocity1[2]
    
    bead2diff = [0,0,0]
    bead2diff[0] = md_coef2[3]*velocity2[0]
    bead2diff[1] = md_coef2[3]*velocity2[1]
    bead2diff[2] = md_coef2[3]*velocity2[2]
    
    bead1 = bead1+bead1diff
    bead2 = bead2+bead2diff
    
print('### Initial Positions ###')
print('Bead1: ',bead1start)
print('Bead2: ',bead2start)    
print('### Final Positions ###')
print('Bead1: ',bead1)
print('Bead2: ',bead2)
print('### Difference ###')
print('bead1: ', bead1-bead1start)
print('bead2: ', bead2-bead2start)

'''
Prompt

This question relates to molecular dynamics using Langevin dynamics with velocities. Two atoms are covalently bonded with a harmonic potential of:
U_{bond} = \frac{1}{2} \sum_i k (d - r0)^2
Where k is a force constant of 15.0 kcal/mol/angstrom^2 and r0 = 6.13 angstrom. The beads start at positions: 
bead1: 743.2363, 45.6731, 153.5673
bead2: 741.0813, 48.6721, 158.6731
where a random bump acts upon them to move them. Bead1 has a bump of [1.2856054152322811, -0.30355337730615073, 0.6190756618694475], and bead2 has a bump of [0.39599854914377247, 0.22340564775539556, -0.05433941722861124]. dT of this simulation is 0.2 and as Langevin dynamics are used, there is not unit. Friction set at 0.0001. The mass of bead 1 is 328.212 Da and bead 2 is 344.212 Da. The system maintained at 298.15 K. Forces are in kcal/mol/angstrom and positions are in angstrom. 
The MD coefficients for integration can with the following function
def set_MD_Coef():
    md_coef1 = [0,0,0,0]
    c1 = 0.5 * dt * friction / mass1
    c2 = np.sqrt(1.0 / (1.0 + c1))
    md_coef1[0] = 0.5 * c2 / mass * np.sqrt(2.0 * friction * 0.0019872042586408316 * temperature * dt)
    md_coef1[1] = (1.0 - c1) / (1.0 + c1)
    md_coef1[2] = c2 * dt / mass
    md_coef1[3] = c2 * dt
return md_coef1
What is the coordinate position after 1 simulation step? Print the answer with 6 values in a list as follows: 
"Initial Positions"
"Bead1": [x,y,z]
Bead2": [x,y,z]
"Final Positions"
"Bead1": [x,y,z]
"Bead2": [x,y,z]
"Difference"
"Bead1": [x,y,z]
"Bead2": [x,y,z]

'''