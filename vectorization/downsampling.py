import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans

import os
from tqdm import tqdm
import pickle
import random
import shutil

cluster_number = 4
fpsize = 1024
target_size = 900
dimension = fpsize
random.seed(10)

current_path = os.getcwd()
os.chdir("..")
repo_path = os.getcwd()
os.chdir(current_path)

with open('./data/rdk_fp_sdf_1024.pkl', 'rb') as f:
     fp = pickle.load(f)
df_full = pd.DataFrame(fp).T
df = df_full.dropna()   


pca = PCA(n_components = dimension) 
decomp = pca.fit_transform(df) 
df_pca= pd.DataFrame(decomp)
df_pca.index = df.index

var=np.cumsum(pca.explained_variance_ratio_*100)
index = np.where(var>90)
threshold = index[0][0]
pc_range = np.arange(threshold)

'''labelling based on clustering/PCA'''
kmeans = KMeans(n_clusters=cluster_number, random_state=10)
clusters = kmeans.fit_predict(df_pca[pc_range])

df_pca['clustered'] = pd.Series(clusters, index=df.index)

for i in tqdm(np.arange(cluster_number)):
    df_sub = df[df_pca['clustered'] == i]
    code = df_sub.index

    size = int(np.round(target_size * len(code)/len(df)))
    print(f"size of subgroup{i+1}: {size}")
    subgroup = random.choices(code, k=size)

    '''copy and paste files'''
    xyz_file = []
    for code in subgroup:
        file = code+".xyz"
        xyz_file.append(file)
    
    # copy candidate xyz files
    src_path = repo_path+"/filtered_xyz_data/"
    dst_path = repo_path+"/subset_xyz_data/"
    for file in xyz_file:
        shutil.copy(src_path+file, dst_path)