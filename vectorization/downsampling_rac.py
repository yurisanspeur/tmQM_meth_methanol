import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from sklearn.preprocessing import MaxAbsScaler

import os
from tqdm import tqdm
import pickle
import random
import shutil

'''define path'''
current_path = os.getcwd()
os.chdir("..")
repo_path = os.getcwd()
os.chdir(current_path)

'''load data'''
with open('./data/autocorrelation_oxo_intermed.pkl', 'rb') as f:
     fp = pickle.load(f)
df_full = pd.DataFrame(fp).T
df = df_full.dropna()

'''normalization'''
scaler = MaxAbsScaler()
scaler.fit(df)
scaled = scaler.transform(df)
df_nrm = pd.DataFrame(scaled, columns=df.columns, index=df.index)

cluster_number = 4
target_size = 900
dimension = len(df_nrm.iloc[0])
random.seed(10)

'''PCA'''
pca = PCA(n_components = dimension) 
decomp = pca.fit_transform(df_nrm) 
df_pca= pd.DataFrame(decomp)
df_pca.index = df_nrm.index

var=np.cumsum(pca.explained_variance_ratio_*100)
index = np.where(var>90)
threshold = index[0][0]
pc_range = np.arange(threshold)

'''labelling based on clustering/PCA'''
kmeans = KMeans(n_clusters=cluster_number, random_state=10)
clusters = kmeans.fit_predict(df_pca[pc_range])

df_pca['clustered'] = pd.Series(clusters, index=df_nrm.index)

for i in np.arange(cluster_number):
    df_sub = df_nrm[df_pca['clustered'] == i]
    code = df_sub.index

    size = int(np.round(target_size * len(code)/len(df_nrm)))
    print(f"size of subgroup{i+1}: {size}")
    subgroup = random.sample(list(code), size)
    #print(f"cluster{i+1}: {len(subgroup)}")

    '''copy and paste files'''
    xyz_file = []
    for code in subgroup:
        file = code+".xyz"
        xyz_file.append(file)
    with open('./data/subgroup'+str(i+1)+'_ids.pkl', 'wb') as f:
        pickle.dump(xyz_file, f)
    
    '''copy candidate xyz files'''
    src_path = repo_path+"/individual_xyz_data/"
    dst_path = repo_path+"/subset_xyz_rac/"
    for file in xyz_file:
        shutil.copy(src_path+file, dst_path)
        
'''final output check'''
files = os.listdir(repo_path+"/subset_xyz_rac/")
print(f"the total number of subset: {len(files)}")