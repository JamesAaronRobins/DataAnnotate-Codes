import numpy as np

# --- Constants ---
k = 15.0  # kcal/mol/angstrom^2
r0 = 6.13  # angstrom
dt = 0.2
friction = 0.0001
temperature = 298.15  # K
R_kcal = 0.0019872042586408316 # Gas constant in kcal/mol/K

# --- Bead Properties ---
mass1 = 328.212  # Da
mass2 = 344.212  # Da

# --- Initial Positions (as numpy arrays for vector math) ---
pos1_initial = np.array([743.2363, 45.6731, 153.5673])
pos2_initial = np.array([741.0813, 48.6721, 158.6731])

# --- Random Bumps (as numpy arrays) ---
bump1 = np.array([1.2856054152322811, -0.30355337730615073, 0.6190756618694475])
bump2 = np.array([0.39599854914377247, 0.22340564775539556, -0.05433941722861124])

# --- Initial Velocities (assumed to be zero) ---
vel1 = np.array([0.0, 0.0, 0.0])
vel2 = np.array([0.0, 0.0, 0.0])

# Calculate the vector and distance between the beads
r_vec = pos1_initial - pos2_initial
d = np.linalg.norm(r_vec)

# Calculate the magnitude of the force
# F = -k * (d - r0)
# A negative magnitude means an attractive (restoring) force
force_magnitude = -k * (d - r0)

# Calculate the unit vector along the bond
unit_vector = r_vec / d

# Distribute the force onto each bead (Newton's 3rd Law)
# Force on bead 1 is in the direction of the unit vector
force1 = force_magnitude * unit_vector
# Force on bead 2 is equal and opposite
force2 = -force1

def set_MD_Coef(mass, dt, friction, temperature):
    """
    Calculates the coefficients for the Langevin dynamics integrator.
    """
    md_coef = [0,0,0,0]
    c1 = 0.5 * dt * friction / mass
    c2 = np.sqrt(1.0 / (1.0 + c1))
    # Coefficient for the random force component
    md_coef[0] = 0.5 * c2 / mass * np.sqrt(2.0 * friction * R_kcal * temperature * dt)
    # Coefficient for velocity decay due to friction
    md_coef[1] = (1.0 - c1) / (1.0 + c1)
    # Coefficient for the systematic force component
    md_coef[2] = c2 * dt / mass
    # Coefficient for the position update
    md_coef[3] = c2 * dt
    return md_coef

# Calculate coefficients for each bead
md_coef1 = set_MD_Coef(mass1, dt, friction, temperature)
md_coef2 = set_MD_Coef(mass2, dt, friction, temperature)

# --- Update Bead 1 ---
# v_new = v_old * C1 + F * C2 + bump * C0
vel1_new = vel1 * md_coef1[1] + force1 * md_coef1[2] + bump1 * md_coef1[0]
# x_new = x_old + v_new * C3
pos1_final = pos1_initial + vel1_new * md_coef1[3]

# --- Update Bead 2 ---
vel2_new = vel2 * md_coef2[1] + force2 * md_coef2[2] + bump2 * md_coef2[0]
pos2_final = pos2_initial + vel2_new * md_coef2[3]

# --- Calculate Differences ---
diff1 = pos1_final - pos1_initial
diff2 = pos2_final - pos2_initial

print('"Initial Positions"')
print(f'"Bead1": {pos1_initial.tolist()}')
print(f'"Bead2": {pos2_initial.tolist()}')
print('"Final Positions"')
print(f'"Bead1": {pos1_final.tolist()}')
print(f'"Bead2": {pos2_final.tolist()}')
print('"Difference"')
print(f'"Bead1": {diff1.tolist()}')
print(f'"Bead2": {diff2.tolist()}')