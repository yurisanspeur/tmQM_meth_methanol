import pickle
import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
import matplotlib.pyplot as plt

## data load
with open('./data/rdk_fp_sdf_1024.pkl', 'rb') as f:
     rdkfp = pickle.load(f)
df_rdkfp_full = pd.DataFrame(rdkfp).T
df_rdkfp_sdf = df_rdkfp_full.dropna()   
df_rdkfp_sdf

with open('./data/morgan_fp_sdf_1024.pkl', 'rb') as f:
     mgfp = pickle.load(f)
df_mgfp_full = pd.DataFrame(mgfp).T
df_mgfp_sdf = df_mgfp_full.dropna()   
df_mgfp_sdf


dimension = 2
cluster_num = 4

'''PCA'''
pca1 = PCA(n_components = dimension) 
decomp1 = pca1.fit_transform(df_rdkfp_sdf)
ex_var1 = np.sum(pca1.explained_variance_ratio_)
print(decomp1)
print(f'RDKit fp explained variance: {ex_var1}')

pca2 = PCA(n_components = dimension) 
decomp2 = pca2.fit_transform(df_mgfp_sdf) 
ex_var2 = np.sum(pca2.explained_variance_ratio_)
print(decomp2)
print(f'Morgan fp explained variance: {ex_var2}')

'''K-Means'''
kmeans1 = KMeans(n_clusters = cluster_num) 
kmeans1.fit(df_rdkfp_sdf) 

kmeans2 = KMeans(n_clusters = cluster_num)  
kmeans2.fit(df_mgfp_sdf) 

x1 = decomp1[:, 0] 
y1 = decomp1[:, 1]

x2 = decomp2[:, 0] 
y2 = decomp2[:, 1]

 
# rdkit
# plt.subplot(1,2,1)
plt.figure(figsize = (6,5))
plt.scatter(x1, y1, c = kmeans1.labels_, alpha = 0.7) 
plt.title("PCA: rdkit_fps, cluster") 
plt.colorbar()
plt.savefig("./plots/visual_rdkit_k"+str(cluster_num)+".png")
# morgan
# plt.subplot(1,2,2) 
plt.figure(figsize = (6,5))
plt.scatter(x2, y2, c = kmeans2.labels_, alpha = 0.7)
plt.title("PCA: morgan_fps, cluster") 
plt.colorbar()
plt.savefig("./plots/visual_morgan_k"+str(cluster_num)+".png")