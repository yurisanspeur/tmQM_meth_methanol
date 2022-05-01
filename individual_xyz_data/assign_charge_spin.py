from pymatgen.core import Molecule
import pandas as pd
import sys
# Mapping for orig charge and multiplicity

mol = sys.argv[1]

df = pd.read_csv('charge_mult_orig.csv')
orig_charge = df[df['CSD_code'] == mol]['charge'].iloc[0]
orig_mnd = df[df['CSD_code'] == mol]['nn'].iloc[0]
orig_spin = 0


import pdb; pdb.set_trace()
mol = Molecule.from_file(f"target_intermediate_{mol}_optimised_xtb.xyz")

mol_comp = mol.composition
target_charge = orig_charge + 1 # in case of methanol
mol_comp_ox = mol_comp.add_charges_from_oxi_state_guesses(target_charge=target_charge) # Get the charge from the csv and account for fact that we pop halide and the charge on the intermediates
# Get the metal site
metal_center = [el for el in mol_comp_ox.elements if el.is_metal][0]

# Get the spin from the coordination (take from MND), the oxidation state and the location in the periodic table
if metal_center.row >= 5:
    spin_config = "low"
else:
    spin_config = "high"

if orig_mnd > 5:
    coordination = 'oct'
else:
    coordination = 'tet'

spin = metal_center.get_crystal_field_spin(coordination=coordination, spin_config=spin_config)

