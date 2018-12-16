import matplotlib.pyplot as plt
import numpy as np
import sys

#Macros
F = 1000
k = 50
epsilon = 0.5

#  randomly initiate k local and global optimums in a 2-D space
# local_optimum = np.random.rand(2,k)*F
# global_optimum = np.random.rand(2,k)*F
# np.savetxt("local_optimum.txt",local_optimum)
# np.savetxt("global_optimum.txt",global_optimum)
local_optimum = np.loadtxt('local_optimum.txt')
global_optimum = np.loadtxt('global_optimum.txt')

# print (global_optimum)
plt.scatter(local_optimum[0],local_optimum[1],facecolors='orange')
plt.scatter(global_optimum[0],global_optimum[1],facecolors='royalblue')
plt.show()

# calculate Di and Distar
# Di is the closest global optimum to a local optimum,
# Distar the closest local optimum to a global optimum
Di = np.zeros(k)
Distar = np.zeros(k)

for i in range(k):
    closest = -1
    d = 1000000
    for j in range(k):
        curr_distance = np.sqrt((global_optimum[0][j]-local_optimum[0][i])**2 + (global_optimum[1][j]-local_optimum[1][i])**2)
        if curr_distance < d:
            d = curr_distance
            closest = j
    Di[i] = d

for i in range(k):
    closest = -1
    d = 1000000
    for j in range(k):
        curr_distance = np.sqrt((global_optimum[0][i]-local_optimum[0][j])**2 + (global_optimum[1][i]-local_optimum[1][j])**2)
        if curr_distance < d:
            d = curr_distance
            closest = j
    Distar[i] = d

# global_bar[i] = j, j is the index the point in Globaloptimum 
Global_bar = []
Local_bar = []
#  we need to have a non-descending of Di and D_i*
Di_ascending_index = np.argsort(Di)
Distar_ascending_index = np.argsort(Distar)

#  global_sparse_mapping[i] = j means that η(i) = j
Global_sparse_mapping = np.zeros(k)
Local_sparse_mapping = np.zeros(k)

for index in range(k):
    i = Distar_ascending_index[index]
    flag = False
    # if ∃ i ∗ ∈ O such that δ(i ∗ , i ∗ bar ) ≤ epsilon · D i ∗ then
    for j in Global_bar:
        if np.sqrt((global_optimum[0][j]-global_optimum[0][i])**2 + (global_optimum[1][j]-global_optimum[1][i])**2) <= epsilon * Distar[i]:
            Global_sparse_mapping[i] = j
            flag = True
            break
    # else
    if not flag:
        Global_sparse_mapping[i] = i
        Global_bar.append(i)

for index in range(k):
    # i is i
    i = Di_ascending_index[index]
    flag = False
    # if ∃ i  ∈ O such that δ(i , i  bar ) ≤ epsilon · D i  then
    #  j is bar i 
    for j in Local_bar:
        if np.sqrt((local_optimum[0][j]-local_optimum[0][i])**2 + (local_optimum[1][j]-local_optimum[1][i])**2) <= epsilon * Di[i]:
            Local_sparse_mapping[i] = j
            flag = True
            break
    # else
    if not flag:
        Local_sparse_mapping[i] = i
        Local_bar.append(i)

plt.scatter(local_optimum[0][Local_bar],local_optimum[1][Local_bar],facecolors='orange',label='local_bar')
plt.scatter(global_optimum[0][Global_bar],global_optimum[1][Global_bar],facecolors='royalblue',label= 'global_bar')

local_filtered_out = np.setdiff1d(np.arange(k),Local_bar)
global_filtered_out = np.setdiff1d(np.arange(k),Global_bar)
plt.scatter(local_optimum[0][local_filtered_out],local_optimum[1][local_filtered_out],facecolors='white',edgecolor='orange',label='local_filtered out')
plt.scatter(global_optimum[0][global_filtered_out],global_optimum[1][global_filtered_out],facecolors='white',edgecolor='royalblue',label='global_filtered_out')
plt.legend()
plt.savefig("Sparseed.png")
plt.show()

# Map the each point in global_bar or local_bar to its nearest neighbor with diff color
# local_phi_bar[i] = j represents local_optimal[local_bar[i]]'s nearset neighbor is global_optimal[global_bar[j]]
Local_closest_bar = np.zeros(len(Local_bar))
Global_closest_bar = np.zeros(len(Global_bar))

# local_center_mapping[i] = B means B is the closest global optimum among all points in global bar whose closest local optimum in local_bar is local_bar[i]
# global_center_mapping[i] = B means B is the closest local optimum among all points in local bar whose closest global optimum in global_bar is global_bar[i]
# they map to the address in local_bar or global_bar instead of the actual index
local_center_mapping = {}
global_center_mapping = {}
# find center(i) for i* in global bar 
for i in range(len(Local_bar)):
    d_min = 1000000
    closest = -1
    for j in range(len(Global_bar)):
        curr_distance = np.sqrt( (local_optimum[0][Local_bar[i]] - global_optimum[0][Global_bar[j]])**2 + (local_optimum[1][Local_bar[i]] - global_optimum[1][Global_bar[j]])**2)
        if curr_distance < d_min:
            d_min = curr_distance
            closest = j
    Local_closest_bar[i] = d_min
    # print (d_min)
    # print (closest)
    # print (local_optimum[0][Local_bar[i]])
    # print (global_optimum[0][Global_bar[closest]])
    # print (local_optimum[1][Local_bar[i]])
    # print (global_optimum[1][Global_bar[closest]])
    # print (np.sqrt( (local_optimum[0][Local_bar[i]] - global_optimum[0][Global_bar[closest]])**2 + (local_optimum[1][Local_bar[i]] - global_optimum[1][Global_bar[closest]])**2))
    if closest not in global_center_mapping:
        global_center_mapping[closest] = [(i,d_min)]
    else:
        global_center_mapping[closest].append((i,d_min))
    # exit()


for k,v in global_center_mapping.items():
    global_center_mapping[k] = min(v,key = lambda x:x[1])
    
#  find center(i) for i in local bar 
for i in range(len(Global_bar)):
    d_min = 1000000
    closest = -1
    for j in range(len(Local_bar)):
        curr_distance = np.sqrt( (local_optimum[0][Local_bar[j]] - global_optimum[0][Global_bar[i]])**2 + (local_optimum[1][Local_bar[j]] - global_optimum[1][Global_bar[i]])**2)
        if curr_distance < d_min:
            d_min = curr_distance
            closest = j
    Global_closest_bar[i] = d_min
    if closest not in local_center_mapping:
        local_center_mapping[closest] = [(i,d_min)]
    else:
        local_center_mapping[closest].append((i,d_min))

for k,v in local_center_mapping.items():
    local_center_mapping[k] = min(v,key = lambda x:x[1])

# global_center_mapping[i] = B means B is the closest local optimum among all points in local bar whose closest global optimum is global_bar[i] 
# now we want to a graph which connect i in and cent(i) so we can see visualization of T

# what we have now: global_center_map[i] = j means that local_bar[j] is the most close one whose center is global_bar[i]
for k,v in local_center_mapping.items():
    pt1 = Local_bar[k]
    pt2 = Global_bar[v[0]]
    plt.plot([local_optimum[0][pt1],global_optimum[0][pt2]],[local_optimum[1][pt1],global_optimum[1][pt2]],color='m')

plt.scatter(local_optimum[0][Local_bar],local_optimum[1][Local_bar],facecolors='orange',label='local_bar')
plt.scatter(global_optimum[0][Global_bar],global_optimum[1][Global_bar],facecolors='royalblue',label= 'global_bar')

plt.scatter(local_optimum[0][local_filtered_out],local_optimum[1][local_filtered_out],facecolors='white',edgecolor='orange',label='local_filtered out')
plt.scatter(global_optimum[0][global_filtered_out],global_optimum[1][global_filtered_out],facecolors='white',edgecolor='royalblue',label='global_filtered_out')
plt.legend()
plt.savefig("Tunfiltered.png")
plt.show()

T = {}
# There is additional requirement for T.
# That T i  s required epsilon · δ(cent(i), i) ≤ D i
for k,v in local_center_mapping.items():
    pt1 = Local_bar[k]
    pt2 = Global_bar[v[0]]
    dist = v[1]
    if dist *epsilon <= Di[pt1]:
        T[k] = v

for k,v in T.items():
    pt1 = Local_bar[k]
    pt2 = Global_bar[v[0]]
    plt.plot([local_optimum[0][pt1],global_optimum[0][pt2]],[local_optimum[1][pt1],global_optimum[1][pt2]],color='m')

plt.scatter(local_optimum[0][Local_bar],local_optimum[1][Local_bar],facecolors='orange',label='local_bar')
plt.scatter(global_optimum[0][Global_bar],global_optimum[1][Global_bar],facecolors='royalblue',label= 'global_bar')

plt.scatter(local_optimum[0][local_filtered_out],local_optimum[1][local_filtered_out],facecolors='white',edgecolor='orange',label='local_filtered out')
plt.scatter(global_optimum[0][global_filtered_out],global_optimum[1][global_filtered_out],facecolors='white',edgecolor='royalblue',label='global_filtered_out')
plt.legend()
plt.savefig("Tfiltered.png")
plt.show()

#  N is the net discussed in section 1.2
N = set()
for local_index in Local_bar:
    for global_index in Global_bar:
        distance = np.sqrt( (global_optimum[0][global_index]-local_optimum[0][local_index])**2 + (global_optimum[1][global_index]-local_optimum[1][local_index])**2     )
        if distance<= Di[local_index]/epsilon and Distar[global_index]>= epsilon* Di[local_index]:
            N.add((global_index,local_index))

for pair in N:
    global_index = pair[0]
    local_index = pair[1]
    plt.plot([local_optimum[0][local_index],global_optimum[0][global_index]],[local_optimum[1][local_index],global_optimum[1][global_index]],color='lawngreen')

plt.scatter(local_optimum[0][Local_bar],local_optimum[1][Local_bar],facecolors='orange',label='local_bar')
plt.scatter(global_optimum[0][Global_bar],global_optimum[1][Global_bar],facecolors='royalblue',label= 'global_bar')

plt.scatter(local_optimum[0][local_filtered_out],local_optimum[1][local_filtered_out],facecolors='white',edgecolor='orange',label='local_filtered out')
plt.scatter(global_optimum[0][global_filtered_out],global_optimum[1][global_filtered_out],facecolors='white',edgecolor='royalblue',label='global_filtered_out')
plt.legend()
plt.savefig("N.png")
plt.show()
