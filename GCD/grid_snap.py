# Ver 1.0 of the grid fitting algorithm
# At this stage, it assumes that the grids are sufficiently far away from each
# other that an optimal position will always be found for each of the clusters.
import numpy as np
import csv
# How does it work?
# This algorithm takes the clusters and, treating each string as an anchor,
# looks for the best possible set of transformations to snap it to grid.

# Rotates the clusters so that the fitting may be tested 
def rotate(cluster,n):
    angle=n*5*np.pi/180
    R=np.array([[np.cos(angle),-np.sin(angle)],[np.sin(angle),np.cos(angle)]])
    centroid = np.mean(cluster, axis=0)
    return [R @ (p - centroid) + centroid for p in cluster]
# Checks if a given point in the generated geometry is within 40 meters of the
# members of the cluster that have already been found.
def clust(temp,point,spacing):
    for i in temp:
        dist=i-point
        if np.isclose(np.sqrt(dist[0]**2 + dist[1]**2),spacing):
            return True
    return False
# Does the heavy lifting of creating the clusters and arranging them into a 2d array
def clustering(pos,nstring,spacing):
    clusters=[]
    idx=[]
    for i in range(len(pos)):
        if i not in idx:
            idx.append(i)
            temp=[pos[i]]
            while len(temp)!=nstring:
                for j in range(len(pos)):
                    if j not in idx:
                        if clust(temp,pos[j],spacing)==True:
                            temp.append(pos[j])
                            idx.append(j)
            clusters.append(temp)
    print("Clustered!")
    return clusters
# Returns the distance for the displacement vector
def distance(displacement):
    returnal=[]
    for i in displacement:
        returnal.append(np.sqrt(i[0]**2+i[1]**2))
    return returnal
# Snaps a given point to the nearest position in the grid and checks the distortion 
# resulting from it.
def snap(grid,point,new_pos,used):
    displacement=grid-point
    dist=distance(displacement)
    idx=np.argmin(dist)
    while idx in used:
        dist[idx]=100000
        idx=np.argmin(dist)
    if np.count_nonzero(dist==dist[idx])!=1:
        distances=[]
        for i in np.where(dist==dist[idx])[0]:
            disp=new_pos-grid[i]
            distances.append(np.sum(distance(disp)))
        idx=np.argmin(distances)
    used.append(idx)
    return idx,displacement[idx]
# Builds the new clusters considering the optimal anchor and translation
def cluster_builder(idxs,grid):
    cluster=[]
    for i in idxs:
        cluster.append(grid[i])
    return cluster
# Moves a cluster and checks the resulting distortion, then saves all of the 
# possible translations for a given anchor in an array
def anchoring(clusters,grid,index,idxs,tolerance,gamma1,gamma2,gamma3):
    new_clusters=[]
    scores=[]
    n=0
    for i in clusters:
        best=1000000
        for m in range(36):
            i=rotate(i,1)
            displacement= grid-i[index]
            dist=distance(displacement)
            found=False
            while found==False and np.min(dist)<=tolerance:
                skip=False
                idx=np.argmin(dist)
                for j in range(len(clusters)):
                    if j!=n and idx in idxs[j]:
                        dist[idx]=100000
                        skip=True
                if skip==True:
                    continue
                new_pos=np.array(i)+displacement[idx]
                alligned=True
                total=[]
                temp=[]
                used=[]
                for j in new_pos:
                    snaps=snap(grid,j,new_pos,used)
                    total.append(snaps[1])
                    temp.append(snaps[0])
                    if not np.any((grid == j).all(axis=1)):
                        alligned=False
                mean=np.mean(total,axis=0)
                total=np.array(total)
                total=total-mean
                if gamma1*np.sum(distance(total))+gamma2*m+gamma3*np.min(dist)<best:
                    best=gamma1*np.sum(distance(total))+gamma2*m+gamma3*np.min(dist)
                    idxs[n]=np.copy(temp)
                if alligned==True:
                    found=True
                else:
                    dist[idx]=100000
        scores.append(best)
        new_clusters.append(cluster_builder(idxs[n],grid))
        n=n+1
    return new_clusters,scores
# Wrapper function for the whole fitting process
def translation(cluster, grid, gamma1,gamma2,gamma3,tolerance):
    used=[]
    new_cluster=[]
    best=np.full(len(cluster),100000)
    for i in range(len(cluster)):
        used.append([])
        new_cluster.append([])
    for i in range(len(cluster[0])):
        temp_clusters,temp_scores=anchoring(cluster,grid,i,used,tolerance,gamma1,gamma2,gamma3)
        for j in range(len(cluster)):
            if temp_scores[j]<best[j]:
                new_cluster[j]=temp_clusters[j]
                best[j]=temp_scores[j]
    return new_cluster
# Wrapper function that does the clustering and fitting
def fit(grid,pos,nstring,gamma1=1,gamma2=10,gamma3=0,tolerance=40,spacing=40):
    clusters=clustering(pos,nstring,spacing)
    final=translation(clusters,grid,gamma1,gamma2,gamma3,tolerance)
    return np.array(final)
# Just a wrapper to make creating the arrays easier
def group(x,y):
    return np.column_stack((x, y))
# Once again, just a wrapper to make creating the arrays easier
def create_grid(grid_path):
    x=[]
    y=[]
    with open(grid_path, mode='r', newline='', encoding='utf-8') as file:
        reader = csv.reader(file)
        for row in reader:
            y.append(float(row[1]))
            x.append(float(row[0]))
    grid=np.column_stack((x, y))
    return grid
# Flattens the output of the fitting process into two arrays for 
# the x and y coordinates
def flatten(clusters):
    clusters=clusters.reshape(-1,2)
    return clusters[:,0],clusters[:,1]
# How to use:
# - First, call the create_grid function with the path to your grid file (OPTIONAL).
# - Then, provide the group function with two arrays for the x and y coordinates of
#   your grid points (OPTIONAL).
# - With these two, call the fit function. It requires 3 arguments:
#       - grid : Your grid to fit the points to 
#       - pos  : The position of your strings
#       - nstring : The number of strings in each of your clusters
#   The function will return you a 3D array with the new string positions
#   organized in clusters for your convenience.
# - In case your particular application needs the flat coordinates, 
#   call the flatten function.  

#   Penalizations (or the use of gammas):
# Different applications require different things. For your convenience,
# this library provides you with three parameters that control how much
# the distortion, rotation of the cluster and distance of translation are
# penalized
#   - gamma1: Controls the penalization for the distortion (aka the sum of
#     the displacement for each snapped point). Default value is 1.
#   - gamma2: Controls the penalization for the rotation. Default value is
#     10. Lower values will cause your clusters to be lined up!
#   - gamma3: Controls the penalization for the translation of the anchor
#     point. Default value is 0.