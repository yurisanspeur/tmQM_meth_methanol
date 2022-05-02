import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler, MaxAbsScaler
from sklearn.metrics import silhouette_samples, silhouette_score

import matplotlib.pyplot as plt
from matplotlib import gridspec
import matplotlib.cm as cm
from tqdm import tqdm
import pickle


'''data load'''
with open('./data/autocorrelation.pkl', 'rb') as f:
     fp = pickle.load(f)
df_full = pd.DataFrame(fp).T
df = df_full.dropna()   

'''normalization'''
scaler = MaxAbsScaler()
scaler.fit(df)
scaled = scaler.transform(df)
df_nrm = pd.DataFrame(scaled, columns=df.columns, index=df.index)


'''PCA'''
dimension = len(df_nrm.iloc[0])

pca = PCA(n_components = dimension) 
decomp = pca.fit_transform(df_nrm) 
df_pca= pd.DataFrame(decomp)
df_pca.index = df_nrm.index


'''variance explained'''
plt.figure(figsize=(8,6))
fig, ax = plt.subplots(figsize=(8,6))

var=np.cumsum(pca.explained_variance_ratio_*100)
index = np.where(var>90)
threshold = index[0][0]

plt.plot([i+1 for i in range(len(var))],var,'-',linewidth=2)
plt.xticks([i+1 for i in range(len(var)) if i==0 or i%100 == 99])
plt.ylabel('% Variance Explained',fontsize=16,fontweight='bold')
plt.xlabel('Pincipal Component (PC)',fontsize=16,fontweight='bold')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.tight_layout()
plt.tick_params ('both',width=2,labelsize=12)
plt.axvline(threshold, ls='--', c='k');
plt.savefig("./plots/rac_var.png",
            bbox_inches='tight', 
             facecolor='w');
print(f'explained variance > 90% at {threshold+1:.1f}')


'''sum of distance'''
sse = []
list_k = list(range(1, 10))
pc_range = np.arange(threshold+1)

for k in tqdm(list_k):
    km = KMeans(n_clusters=k)
    km.fit(df_pca[pc_range])
    sse.append(km.inertia_)

# Plot sse against k
plt.figure(figsize=(8,6))
plt.plot(list_k, sse, '-o')
plt.xlabel(r'Number of clusters k', fontsize=16,fontweight='bold')
plt.ylabel('Sum of squared distance', fontsize=16,fontweight='bold')
plt.savefig("./plots/rac_sum_sqdist.png",
            bbox_inches='tight', 
             facecolor='w');



'''silhouette'''
range_n_clusters = [2, 3, 4, 5, 6, 7, 8, 9, 10]
for n_clusters in range_n_clusters:
    fig, (ax1,ax2)= plt.subplots(1, 2)
    fig.set_size_inches(8, 4)
    
    kmeans = KMeans(n_clusters=n_clusters, random_state=10)
    cluster_labels = kmeans.fit_predict(df_pca[pc_range])
    silhouette_avg = silhouette_score(df_pca[pc_range], cluster_labels)
    print("For n_clusters =", n_clusters,
          "The average silhouette_score is :", silhouette_avg)
    
    sample_silhouette_values = silhouette_samples(df_pca[pc_range], cluster_labels)

    y_lower = 10
    
    for i in range(n_clusters):
        # Aggregate the silhouette scores for samples belonging to
        # cluster i, and sort them
        ith_cluster_silhouette_values = \
            sample_silhouette_values[cluster_labels == i]

        ith_cluster_silhouette_values.sort()

        size_cluster_i = ith_cluster_silhouette_values.shape[0]
        y_upper = y_lower + size_cluster_i

        color = cm.nipy_spectral(float(i) / n_clusters)
        ax1.fill_betweenx(np.arange(y_lower, y_upper),
                          0, ith_cluster_silhouette_values,
                          facecolor=color, edgecolor=color, alpha=0.7)

        # Label the silhouette plots with their cluster numbers at the middle
        ax1.text(-0.05, y_lower + 0.5 * size_cluster_i, str(i))

        # Compute the new y_lower for next plot
        y_lower = y_upper + 10  # 10 for the 0 samples

    ax1.set_title("The silhouette plot for the various clusters.")
    ax1.set_xlabel("The silhouette coefficient values")
    ax1.set_ylabel("Cluster label")

    # The vertical line for average silhouette score of all the values
    ax1.axvline(x=silhouette_avg, color="red", linestyle="--")


    # 2nd Plot showing the actual clusters formed
    colors = cm.nipy_spectral(cluster_labels.astype(float) / n_clusters)
    ax2.scatter(df_pca[0], df_pca[1], 
                marker='.', s=30, lw=0, alpha=0.7,c=colors, edgecolor='k')


    # Labeling the clusters
    centers = kmeans.cluster_centers_
    # Draw white circles at cluster centers
    ax2.scatter(centers[:, 0], centers[:, 1], marker='o',
                c="white", alpha=1, s=200, edgecolor='k')

    for i, c in enumerate(centers):
        ax2.scatter(c[0], c[1], marker='$%d$' % i, alpha=1,
                    s=50, edgecolor='k')

    ax2.set_title("The visualization of the clustered data.")
    ax2.set_xlabel("PC1")
    ax2.set_ylabel("PC2")

    plt.suptitle(("Silhouette analysis for KMeans clustering on sample data "
                  "with n_clusters = %d" % n_clusters),
                 fontsize=14, fontweight='bold')
    plt.savefig("./plots/rac_silhouette_var90_k"+str(n_clusters)+".png",
               bbox_inches='tight', 
               facecolor='w')
    
#plt.show()