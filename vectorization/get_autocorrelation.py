import os
import numpy as np
import warnings
from tqdm.auto import tqdm
import pickle
import pandas as pd

from molSimplify.Informatics.autocorrelation import *
from molSimplify.Informatics.misc_descriptors import *
from molSimplify.Informatics.graph_analyze import *

def make_rac(xyz_file, m_depth, l_depth, is_oct):
    properties = ['electronegativity', 'nuclear_charge', 'ident', 'topology', 'size', 
                'group_number', 'polarizability', 'num_bonds'] 
    this_mol = mol3D()
    this_mol.readfromxyz(xyz_file)
    feature_names = []
    mc_corrs = np.zeros(shape=(len(properties), (m_depth+1)))
    metal_idx = this_mol.findMetal()[0]
    mc_delta_metricz =  np.zeros(shape=(len(properties), m_depth))
    for idx, p in enumerate(properties):
        delta_list = list(np.asarray(atom_only_deltametric(this_mol, p, m_depth, metal_idx, oct=is_oct)).flatten())
        del delta_list[0]
        mc_corrs[idx] = np.asarray(atom_only_autocorrelation(this_mol, p, m_depth, metal_idx, oct=is_oct)).flatten()
        mc_delta_metricz[idx] = delta_list
        
    if is_oct:
        num_connectors = 6
    else:
        num_connectors = 5
    distances = []
    origin = this_mol.coordsvect()[metal_idx]
    for xyz in this_mol.coordsvect():
        distances.append(np.sqrt((xyz[0]-origin[0])**2+(xyz[1]-origin[1])**2+(xyz[2]-origin[2])**2))

    nearest_neighbours = np.argpartition(distances, num_connectors)
    nn = [x for x in nearest_neighbours[:num_connectors+1] if x != 0] 
    rest_of_autoz = np.zeros(shape=(len(properties), l_depth+1))
    rest_of_deltas = np.zeros(shape=(len(properties), l_depth))
    for idx, p in enumerate(properties):
        rest_of_autoz[idx] = atom_only_autocorrelation(this_mol, p, l_depth, nn, oct=is_oct)
        rest_of_deltas[idx] = atom_only_deltametric(this_mol, p, l_depth, nn)[1:]
    rac_res = np.concatenate((mc_corrs, mc_delta_metricz, rest_of_autoz, rest_of_deltas), axis=None)
    return rac_res

if __name__ == "__main__":
    current_path = os.getcwd()
    os.chdir("..")
    repo_path = os.getcwd()
    os.chdir(current_path)
    
    files = os.listdir(repo_path+'/converted_sdf')
    code_list = []
    for file in files:
        code = file.split('.')[0]
        code_list.append(code)
    '''
    variable for rac
    '''
    m = 5
    l = 2
    
    rac_results = {}
    for code in tqdm(code_list):
        file = repo_path+'/filtered_xyz_data/'+code+".xyz"
        try:
            rac = make_rac(xyz_file=file, m_depth=m, l_depth=l, is_oct=True)
        except:
            rac = None
        rac_results.update({code: rac})
        

    #'./data/autocorrelation_m'+str(m)+'_l'+str(l)+'.pkl'
    with open('./data/autocorrelation.pkl', 'wb') as f:
         pickle.dump(rac_results, f)
    
    
    
    df_rac_full = pd.DataFrame(rac_results).T
    df_rac = df_rac_full.dropna()
    
    idx_list = list(df_rac_full.index)
    for i in df_rac.index:
        idx_list.remove(i)
        failed_index = idx_list

    print('autocorrelation conversion results')
    print(f'total data number: {len(df_rac_full)}')
    print(f'failed number: {len(failed_index)}')