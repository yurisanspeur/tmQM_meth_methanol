from openbabel import pybel
from rdkit import Chem
from rdkit import RDPaths
from rdkit.Chem import AllChem
from rdkit.Chem import DataStructs
from rdkit.Chem.PandasTools import LoadSDF
from tqdm import tqdm

import os
import numpy as np
import warnings
warnings.filterwarnings("ignore")

def xyz_to_smiles(fname: str) -> str:
    mol = next(pybel.readfile("xyz", fname))
    smi = mol.write(format="smi")
    return smi.split()[0].strip()

def read_xyz_file(filename, look_for_charge=True):
    atomic_symbols = []
    xyz_coordinates = []
    charge = 0
    title = ""

    with open(filename, "r") as file:
        for line_number, line in enumerate(file):
            if line_number == 0:
                num_atoms = int(line)
            elif line_number == 1:
                title = line
                if "q=" in line:
                    charge = int(line.split("=")[1])
            else:
                line_data = line.split()
                if len(line_data) != 4:
                    continue
                
                atomic_symbol, x, y, z = line_data
                atomic_symbols.append(atomic_symbol)
                xyz_coordinates.append([float(x), float(y), float(z)])

    atoms = [atom for atom in atomic_symbols]
    return atoms, charge, xyz_coordinates

def xyz_to_fp(code, fpsize, fptype):
    current_path = os.getcwd()
    os.chdir("..")
    repo_path = os.getcwd()
    os.chdir(current_path)
    path = repo_path+'/filtered_xyz_data/'

    if fptype == "rdkit":
        try: 
            s = xyz_to_smiles(path+code+".xyz")
            m = Chem.MolFromSmiles(s)
            fp = Chem.RDKFingerprint(m, fpSize = fpsize)
        except:
            fp = None
    if fptype == "morgan":
        try: 
            s = xyz_to_smiles(path+code+".xyz")
            m = Chem.MolFromSmiles(s)
            fp = AllChem.GetMorganFingerprintAsBitVect(m, 2, fpsize)
        except:
            fp = None
    return np.array(fp)

def sdf_to_fp(code, fpsize, fptype):
    current_path = os.getcwd()
    os.chdir("..")
    repo_path = os.getcwd()
    os.chdir(current_path)
    path = repo_path+'/converted_sdf/'
    
    if fptype == "rdkit":
        try: 
            df = LoadSDF(path+code+'.sdf', smilesName='SMILES')
            s = df['SMILES'].to_numpy()[0]
            m = Chem.MolFromSmiles(s)
            fp = Chem.RDKFingerprint(m, fpSize = fpsize)
        except:
            fp = None
            #pass
    if fptype == "morgan":
        try: 
            df = LoadSDF(path+code+'.sdf', smilesName='SMILES')
            s = df['SMILES'].to_numpy()[0]
            m = Chem.MolFromSmiles(s)
            fp = AllChem.GetMorganFingerprintAsBitVect(m, 2, fpsize)
        except:
            fp = None
    return np.array(fp)

# def get_fp_xyz(path, code):
#     smiles = xyz_to_smiles(path+code+".xyz")
#     m = Chem.MolFromSmiles(smiles)
#     info = {}
#     rdkbi = {}
#     try:
#         mgfp = AllChem.GetMorganFingerprint(m,2,bitInfo=info)
#         rdkfp = Chem.RDKFingerprint(m, maxPath=5, bitInfo=rdkbi)
#         # fp = np.zeros((len(rdkfp))
#         # for i, f in enumerate(rdkfp):
#         #     fp[i] = f
#         # fp = []
#         # for f in rdkfp:
#         #     fp.append(f)
#             #fp = np.array(fp)
#     except:
#         mgfp = None
#         info = None
#         rdkfp = None
#     return mgfp, info, rdkfp

# def get_fp_sdf(path, code):
#     try:
#         df = LoadSDF(path+code+'.sdf', smilesName='SMILES')
#         smiles = df['SMILES'].to_numpy()[0]
#         m = Chem.MolFromSmiles(smiles)
#         info = {}
#         rdkbi = {}
#         mgfp = AllChem.GetMorganFingerprint(m,2,bitInfo=info)
#         rdkfp = Chem.RDKFingerprint(m, maxPath=5, bitInfo=rdkbi)
#         # fp = []
#         # for f in rdkfp:
#         #     fp.append(f)
#     except:
#         mgfp = None
#         info = None
#         rdkfp = None
#     return mgfp, info, rdkfp

    