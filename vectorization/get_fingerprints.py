import os
import numpy as np
import warnings
from tqdm import tqdm
import pandas as pd
import pickle
from analysis.xyz_converter import *

from rdkit import RDLogger
lg = RDLogger.logger()
lg.setLevel(RDLogger.CRITICAL)

import warnings
warnings.filterwarnings("ignore")


current_path = os.getcwd()
os.chdir("..")
repo_path = os.getcwd()
os.chdir(current_path)

files = os.listdir(repo_path+'/converted_sdf')
code_list = []
for file in files:
    #file_list.append(file)
    code = file.split('.')[0]
    code_list.append(code)

fpsize = 1024


'''RDKit fingerprints'''
results = {}
for code in tqdm(code_list):
    results.update({code:sdf_to_fp(code, fpsize, "rdkit")})
with open('./data/rdk_fp_sdf_'+str(fpsize)+'.pkl', 'wb') as f:
     pickle.dump(results, f)
        
'''Morgan fingerprints'''
results = {}
for code in tqdm(code_list):
    results.update({code:sdf_to_fp(code, fpsize, "morgan")})
with open('./data/morgan_fp_sdf_'+str(fpsize)+'.pkl', 'wb') as f:
     pickle.dump(results, f)