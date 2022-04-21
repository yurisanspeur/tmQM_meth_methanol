import autode as ade
import numpy as np
from copy import deepcopy
from autode.input_output import atoms_to_xyz_file
from molfunc import print_combined_molecule
import sys


def rotation_matrix_from_vectors(vec1, vec2):
    """Find the rotation matrix that aligns vec1 to vec2
    :param vec1: A 3d "source" vector
    :param vec2: A 3d "destination" vector
    :return mat: A transform matrix (3x3) which when applied to vec1, aligns it with vec2.
    """
    a, b = (vec1 / np.linalg.norm(vec1)).reshape(3), (
        vec2 / np.linalg.norm(vec2)
    ).reshape(3)
    v = np.cross(a, b)
    c = np.dot(a, b)
    s = np.linalg.norm(v)
    kmat = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
    rotation_matrix = np.eye(3) + kmat + kmat.dot(kmat) * ((1 - c) / (s ** 2))
    return rotation_matrix


# Load the methanol
methanol_orig = ade.Molecule("meth_opt.xyz")
methanol_final = deepcopy(methanol_orig)
# Load the complex
mol = sys.argv[1]
charge = int(sys.argv[2])
complex_ = ade.Molecule(f"{mol}", charge=charge)

# Center the complex with the metal center at the origin

complex_.coordinates = complex_.coordinates - complex_.atoms[0].coord


# Get the location of the halogen that we want to pop out #FIXME Make all the indexing dynamic and robust
halide_atoms_indices = [
    i
    for (i, atom) in enumerate(complex_.atoms)
    if atom.label in ["F", "Cl", "Br", "I", "At"]
]
if len(halide_atoms_indices) > 1:
    pop_index = np.random.choice(halide_atoms_indices)  # Nick's work should provide
else:
    pop_index = halide_atoms_indices[0]
halide_location = complex_.atoms.pop(
    pop_index
).coord  # TODO: Get the index from Nick's work
# Place the oxygen in the intermediate at that location to start

complex_final = deepcopy(complex_)
methanol_final.atoms[0].coord = halide_location
# Control the distance between the oxygen and the metal center

vec = methanol_final.atoms[0].coord - complex_.atoms[0].coord
M_O_distance = np.linalg.norm(vec)
# Let's target an M-O distance of 2.2 Ang
distance_factor = 2.5 / M_O_distance
target_location = complex_.atoms[0].coord + distance_factor * vec

methanol_final.atoms[0].coord = target_location

# Now we need to translate the rest of the atoms in the methanol to be in the coordinate frame of reference
# of the complex
# trans = methanol_final.atoms[0].coord - methanol_orig.atoms[0].coord
# rest_new_coords = np.array(methanol_final.coordinates)
# rest_new_coords[1:,:] = rest_new_coords[1:,:] + trans # shift the remaining methanol_final coords
# methanol_final.coordinates = rest_new_coords


# Add the oxygen to the complex object
complex_.atoms.extend([methanol_final.atoms[0]])

bondlength = 1.429
y_pos = methanol_final.atoms[0].coord[1]
angle = 30  # Initial guess; tune to meet constraints
# Now add the hydrogen such that M-O-H is 120 and O-H bond is 1 Angs
methanol_final.atoms[2].coord = [
    methanol_final.atoms[0].coord[0] + bondlength * np.sin(angle),
    y_pos,
    methanol_final.atoms[0].coord[2] + bondlength * np.cos(angle),
]

complex_.atoms.extend([methanol_final.atoms[2]])


O_index = len(complex_.atoms) - 2
# H_index = len(complex_.atoms) - 1
C_index = len(complex_.atoms) - 1

target_angle_COH = 116

COH_angle = complex_.angle(0, O_index, C_index).to(units="deg")
while COH_angle > target_angle_COH and abs(COH_angle - target_angle_COH) > 0.001:
    angle -= 0.001
    print(angle, "angle")
    methanol_final.atoms[2].coord = [
        methanol_final.atoms[0].coord[0] + bondlength * np.sin(angle),
        y_pos,
        methanol_final.atoms[0].coord[2] + bondlength * np.cos(angle),
    ]
    COH_angle = complex_.angle(0, O_index, C_index).to(units="deg")
    print(COH_angle, "COH_angle")

while COH_angle < target_angle_COH and abs(COH_angle - target_angle_COH) > 0.001:
    angle += 0.001
    print(angle, "angle")
    methanol_final.atoms[2].coord = [
        methanol_final.atoms[0].coord[0] + bondlength * np.sin(angle),
        y_pos,
        methanol_final.atoms[0].coord[2] + bondlength * np.cos(angle),
    ]
    COH_angle = complex_.angle(0, O_index, C_index).to(units="deg")
    print(COH_angle, "COH_angle")


# Now we add the carbon atom with the constraint that the C-O-H angle is 106 degrees and the O-C bond is 1.40 Ang
angle_H = 30
y_pos_H = methanol_final.atoms[1].coord[1]
H_O_bondlength = 1
methanol_final.atoms[1].coord = [
    methanol_final.atoms[0].coord[0] + H_O_bondlength * np.sin(angle_H),
    methanol_final.atoms[0].coord[1],
    methanol_final.atoms[0].coord[2] + H_O_bondlength * np.cos(angle_H),
]  # Hydrogen z same Oxo pivot for x and y

complex_.atoms.extend([methanol_final.atoms[1]])  # Add the H atom

H_index = len(complex_.atoms) - 1
# C_index -= 1
# O_index = len(complex_.atoms) - 3
# O_index -= 1
MOH_angle = complex_.angle(0, O_index, H_index).to(units="deg")
COH_angle = complex_.angle(C_index, O_index, H_index).to(units="deg")
target_angle_MOH = 120
target_angle_COH = 106
# while MOH_angle > target_angle_MOH and abs(MOH_angle - target_angle_MOH) > 0.001:
#    angle_H -= 0.001
#    print(angle_H, "angle_H")
#    methanol_final.atoms[1].coord = [methanol_final.atoms[0].coord[0] + H_O_bondlength*np.sin(angle_H), methanol_final.atoms[0].coord[1], methanol_final.atoms[0].coord[2] + H_O_bondlength*np.cos(angle_H)] #Hydrogen z same Oxo pivot for x and y
#    MOH_angle = complex_.angle(0, O_index, H_index).to(units="deg")
#    print(MOH_angle, "MOH_angle")
#
#
# while MOH_angle < target_angle_MOH and abs(MOH_angle - target_angle_MOH) > 0.001:
#    angle_H += 0.001
#    print(angle_H, "angle_H")
#    methanol_final.atoms[1].coord = [methanol_final.atoms[0].coord[0] + H_O_bondlength*np.sin(angle_H), methanol_final.atoms[0].coord[1], methanol_final.atoms[0].coord[2] + H_O_bondlength*np.cos(angle_H)] #Hydrogen z same Oxo pivot for x and y
#    MOH_angle = complex_.angle(0, O_index, H_index).to(units="deg")
#    print(MOH_angle, "MOH_angle")

while COH_angle > target_angle_COH and abs(COH_angle - target_angle_COH) > 0.001:
    angle_H -= 0.001
    print(angle_H, "angle_H")
    methanol_final.atoms[1].coord = [
        methanol_final.atoms[0].coord[0] + H_O_bondlength * np.sin(angle_H),
        methanol_final.atoms[0].coord[1],
        methanol_final.atoms[0].coord[2] + H_O_bondlength * np.cos(angle_H),
    ]  # Hydrogen z same Oxo pivot for x and y
    COH_angle = complex_.angle(C_index, O_index, H_index).to(units="deg")
    print(COH_angle, "COH_angle")


while COH_angle < target_angle_COH and abs(COH_angle - target_angle_COH) > 0.001:
    angle_H += 0.001
    print(angle_H, "angle_H")
    methanol_final.atoms[1].coord = [
        methanol_final.atoms[0].coord[0] + H_O_bondlength * np.sin(angle_H),
        methanol_final.atoms[0].coord[1],
        methanol_final.atoms[0].coord[2] + H_O_bondlength * np.cos(angle_H),
    ]  # Hydrogen z same Oxo pivot for x and y
    COH_angle = complex_.angle(C_index, O_index, H_index).to(units="deg")
    print(COH_angle, "COH_angle")

# Need a monovalent atom of the same size as carbon

# complex_.atoms[-1].label = 'Be'

# Export the O-H-C backbone to an .xyz file

atoms_to_xyz_file(complex_.atoms, "backbone.xyz")


# Finally, we map the methyl group based on the C atom

# angle_H = -31
# C_methyl_bondlength = 1.1
## Let's try to place a H atom at distance of 1.1 Ang from the C
# methanol_final.atoms[3].coord = [methanol_final.atoms[2].coord[0], methanol_final.atoms[2].coord[1] + C_methyl_bondlength * np.sin(angle_H), methanol_final.atoms[2].coord[2] + C_methyl_bondlength * np.cos(angle_H)]
atoms_to_del = [len(complex_.atoms) - 1]
print_combined_molecule(
    core_xyz_filename="backbone.xyz",
    atoms_to_del=atoms_to_del,
    frag_names=["Me"],
    name="target_intermediate",
)

# rot_matrix = -1 * rotation_matrix_from_vectors(np.array(methanol_orig.atoms[1].coord), np.array(complex_.atoms[H_index].coord))
## Apply rot to methyl atoms
#
# new_methyl_coords = (rot_matrix @ np.array(methanol_orig.coordinates)[3:,:].T).T
#
#
#
#
#
#
##trans = complex_.atoms[C_index].coord - methanol_orig.atoms[2].coord
# rest_new_coords = np.array(methanol_final.coordinates)
# rest_new_coords[3:,:] = new_methyl_coords # map the remaining methanol_final coords
# methanol_final.coordinates = rest_new_coords
#
# complex_final.atoms.extend(methanol_orig.atoms)
# complex_final.atoms.extend(methanol_final.atoms)
# complex_.atoms.extend([methanol_final.atoms[3]])


# from autode.smiles.angles import SDihedral
#
# dihedral = SDihedral(idxs=[26,25,23,24], rot_idxs=[1,0,0,0])
#
# d_angle = dihedral.value(complex_.atoms)
#
# print(d_angle * (180/np.pi))


# C_index = len(complex_.atoms) - 3


# Need a check to see if the M - O - H bond is too small. Let's target 120 degrees
# rot_vec = complex_.atoms[_index].coord - complex_.atoms[O_index].coord # axis of rotation
# rotation_angle = 1
# rot_vec = rot_vec / np.linalg.norm(rot_vec)
#
# while MOH_angle < 120 and abs(MOH_angle - 120) > 0.02:
#    #rotation_angle = 120 - MOH_angle if MOH_angle < 120 else MOH_angle
#    # Unit vector
#    complex_.atoms[H_index].rotate(axis=rot_vec, theta=rotation_angle, origin=complex_.atoms[O_index].coord)
#    MOH_angle = complex_.angle(0, O_index, H_index).to(units="deg")
#    print(MOH_angle)
#
# breakpoint()
#
## Save to xyz file
# atoms_to_xyz_file(complex_final.atoms, "complex_methanol_intermediate.xyz")
# atoms_to_xyz_file(complex_.atoms, "Add_H.xyz")
