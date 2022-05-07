import glob
from tqdm import tqdm
import autode as ade
import numpy as np
from autode.input_output import atoms_to_xyz_file


molecules = glob.glob('*optimised*.xyz')
target_bondlength = 1.84

for mol in tqdm(molecules):
    mol_obj = ade.Molecule(mol, mult=2)
    M_index = np.where(np.array([x.is_metal for x in mol_obj.atoms]))[0][0] 
    curr_MO_distance = mol_obj.distance(M_index, mol_obj.n_atoms - 6)
    factor = target_bondlength / curr_MO_distance
    vec_MO = mol_obj.atoms[mol_obj.n_atoms - 6].coord - mol_obj.atoms[M_index].coord
    target_oxo_pos = mol_obj.atoms[M_index].coord + factor * vec_MO
    translation_vec = target_oxo_pos - mol_obj.atoms[mol_obj.n_atoms - 6].coord
    # Translate the methanol molecule (Note this is a shallow copy so mutating the subset will also mutate it in the complex obj, which we exploit
    methanol = mol_obj.atoms[-6:]
    for atom in methanol:
        atom.translate(translation_vec)

    atoms_to_xyz_file(mol_obj.atoms, f"{mol.split('_')[2]}_translated_methanol_int.xyz")


