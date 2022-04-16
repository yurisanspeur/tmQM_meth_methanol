import autode as ade
import numpy as np
from copy import deepcopy
from autode.input_output import atoms_to_xyz_file

# Load the methanol
methanol_orig = ade.Molecule('methanol_cys.xyz')
methanol_final = deepcopy(methanol_orig)
# Load the complex
complex_ = ade.Molecule('ACUXAX.xyz')

# Center the complex with the metal center at the origin

complex_.coordinates = complex_.coordinates - complex_.atoms[0].coord


# Get the location of the halogen that we want to pop out #FIXME Make all the indexing dynamic and robust
halide_location = complex_.atoms.pop(1).coord
# Place the oxygen in the intermediate at that location to start

methanol_final.atoms[0].coord = halide_location
# Control the distance between the oxygen and the metal center

vec = methanol_final.atoms[0].coord - complex_.atoms[0].coord
M_O_distance = np.linalg.norm(vec)
# Let's target an M-O distance of 2.2 Ang
distance_factor = 1.84 / M_O_distance
target_location = complex_.atoms[0].coord + distance_factor * vec


methanol_final.atoms[0].coord = target_location

# Now we need to translate the rest of the atoms in the methanol to be in the coordinate frame of reference
# of the complex
#trans = methanol_final.atoms[0].coord - methanol_orig.atoms[0].coord
#rest_new_coords = np.array(methanol_final.coordinates)
#rest_new_coords[1:,:] = rest_new_coords[1:,:] + trans # shift the remaining methanol_final coords
#methanol_final.coordinates = rest_new_coords


# Add the oxygen to the complex object
complex_.atoms.extend([methanol_final.atoms[0]]) 

y_pos = methanol_final.atoms[0].coord[1]
angle = 30 # Initial guess; tune to meet constraints
# Now add the hydrogen such that M-O-H is 120 and O-H bond is 1 Angs
methanol_final.atoms[1].coord = [methanol_final.atoms[0].coord[0] + np.sin(angle), y_pos , methanol_final.atoms[0].coord[2] + np.cos(angle)]

complex_.atoms.extend([methanol_final.atoms[1]]) 


O_index = len(complex_.atoms) - 2
H_index = len(complex_.atoms) - 1

MOH_angle = complex_.angle(0, O_index, H_index).to(units="deg")
while MOH_angle > 120 and abs(MOH_angle - 120) > 0.001:
    angle -= 0.001
#    print(angle, "angle")
    methanol_final.atoms[1].coord = [methanol_final.atoms[0].coord[0] + np.sin(angle), y_pos , methanol_final.atoms[0].coord[2] + np.cos(angle)]
    MOH_angle = complex_.angle(0, O_index, H_index).to(units="deg")
#    print(MOH_angle, "MOH_angle")

while MOH_angle < 120 and abs(MOH_angle - 120) > 0.001:
    angle += 0.001
#    print(angle, "angle")
    methanol_final.atoms[1].coord = [methanol_final.atoms[0].coord[0] + np.sin(angle), y_pos , methanol_final.atoms[0].coord[2] + np.cos(angle)]
    MOH_angle = complex_.angle(0, O_index, H_index).to(units="deg")
#    print(MOH_angle, "MOH_angle")

# Now we add the carbon atom with the constraint that the C-O-H angle is 106 degrees and the O-C bond is 1.40 Ang
angle_C = 21
y_pos_H = methanol_final.atoms[1].coord[1]
methanol_final.atoms[2].coord = [methanol_final.atoms[0].coord[0] + 1.4*np.sin(angle_C), methanol_final.atoms[1].coord[1], methanol_final.atoms[0].coord[2] + 1.4*np.cos(angle_C)] #Hydrogen z same Oxo pivot for x and y
breakpoint()
complex_.atoms.extend([methanol_final.atoms[2]])
C_index = len(complex_.atoms) - 1
COH_angle = complex_.angle(C_index, O_index, H_index).to(units="deg")
while COH_angle > 106 and abs(COH_angle - 106) > 0.001:
    angle_C -= 0.001
    print(angle_C, "angle_C")
    methanol_final.atoms[2].coord = [methanol_final.atoms[0].coord[0] + 1.4*np.sin(angle_C), methanol_final.atoms[1].coord[1], methanol_final.atoms[0].coord[2] + 1.4*np.cos(angle_C)] #Hydrogen z same Oxo pivot for x and y
    COH_angle = complex_.angle(C_index, O_index, H_index).to(units="deg")
    print(COH_angle, "COH_angle")


while COH_angle < 106 and abs(COH_angle - 106) > 0.001:
    angle_C += 0.001
    print(angle_C, "angle_C")
    methanol_final.atoms[2].coord = [methanol_final.atoms[0].coord[0] + 1.4*np.sin(angle_C), methanol_final.atoms[1].coord[1], methanol_final.atoms[0].coord[2] + 1.4*np.cos(angle_C)] #Hydrogen z same Oxo pivot for x and y
    COH_angle = complex_.angle(C_index, O_index, H_index).to(units="deg")
    print(COH_angle, "COH_angle")






breakpoint()

#C_index = len(complex_.atoms) - 3


# Need a check to see if the M - O - H bond is too small. Let's target 120 degrees
#rot_vec = complex_.atoms[_index].coord - complex_.atoms[O_index].coord # axis of rotation
#rotation_angle = 1
#rot_vec = rot_vec / np.linalg.norm(rot_vec)
#
#while MOH_angle < 120 and abs(MOH_angle - 120) > 0.02:
#    #rotation_angle = 120 - MOH_angle if MOH_angle < 120 else MOH_angle
#    # Unit vector
#    complex_.atoms[H_index].rotate(axis=rot_vec, theta=rotation_angle, origin=complex_.atoms[O_index].coord)
#    MOH_angle = complex_.angle(0, O_index, H_index).to(units="deg")
#    print(MOH_angle)
#
#breakpoint()
#
## Save to xyz file
atoms_to_xyz_file(complex_.atoms, "complex_methanol_intermediate.xyz")















