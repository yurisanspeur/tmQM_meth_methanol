from pymatgen.core import Molecule
import pandas as pd
import sys
import glob
from tqdm import tqdm
# Mapping for orig charge and multiplicity

molecules = glob.glob('*.xyz')

#mol = sys.argv[1]
for mol in tqdm(molecules):

    try:

        #import pdb; pdb.set_trace()
        mol_name = mol.split('.xyz')[0]

        df = pd.read_csv('charge_mult_oxo.csv')
        orig_charge = df[df['CSD_code'] == mol_name]['charge'].iloc[0]
        orig_mnd = df[df['CSD_code'] == mol_name]['nn'].iloc[0]
        orig_spin = 0


        #import pdb; pdb.set_trace()
        #mol = Molecule.from_file(f"target_intermediate_{mol}_optimised_xtb.xyz")
        mol = Molecule.from_file(f"{mol_name}.xyz_oxo_intermediate.xyz")

        mol_comp = mol.composition
        target_charge = orig_charge # Assume the total charge is the same as the original charge
        #target_charge = orig_charge# in case of oxo

        mol_comp_ox = mol_comp.add_charges_from_oxi_state_guesses(target_charge=target_charge) # The oxidation states are determined from the resting state and then changed based on the
        #intermediate #FIXME: Is this correct?
        # Get the metal site
        metal_center = [el for el in mol_comp_ox.elements if el.is_metal][0]

        # If it is oxo, then we increase the oxidation state by two
        #metal_center._oxi_state += 2 #FIXME: Is this correct?

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
        print(f"Spin for {mol_name} is: {spin+1}")
        # Place charge and mult in df in order to create the orca input files
        df.loc[df['CSD_code'] == mol_name,'oxo_charge'] = int(orig_charge)
        df.loc[df['CSD_code'] == mol_name,'oxo_mult'] = int(spin + 1)
        df.to_csv('charge_mult_oxo.csv', index=False)
    except AttributeError as e:
        print(f"Could not calculate spin for {mol_name}. Failed with exception {e}!")

