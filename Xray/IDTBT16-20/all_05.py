#!/usr/bin/env python
# coding: utf-8

# In[10]:


###########
# Modules #
###########
import numpy as np
import sys, math, os
import MDAnalysis


Atomic_num = dict([ ('C', 6), ('O', 8), ('H', 1), ('S', 16), ('N', 7) ])

WORKDIR = os.getcwd()

snapshot=[1]


binn=0.05

for index, snap in enumerate(snapshot):
    TPRFILE = WORKDIR + "/eq-7.tpr"
    GROFILE = WORKDIR + "/nopbc_"+str(snap)+".gro"

    

    ###########
    #   I/O   #
    ###########
    # Load the universe. In GMX: need for the GRO for the coordinates and the TPR for the connectivity, i.e., topology
    u2 = MDAnalysis.Universe(TPRFILE, GROFILE)
    ALL = u2.select_atoms('not type H')
    
    #--------------------------------------------------------#
    # Extract arrays with atom positions of COM of BT or IDT #
    #--------------------------------------------------------#

    # 0) Make the molecules whole (in place)
    for fragment in ALL.fragments:
        MDAnalysis.lib.mdamath.make_whole(fragment)

    
    ALL_1 = []
    
    for residue in ALL.residues.atoms.positions:
        ALL_1.append(residue)
        
    ALL_1 = np.row_stack(ALL_1).astype('float32')
    ALL_atomicnum = []
    for residue in ALL.residues.atoms.names:
        ALL_atomicnum.append(Atomic_num[residue[0]])
    ALL_1_atomicnum=np.zeros((len(ALL_atomicnum),len(ALL_atomicnum)))
    for i in range(len(ALL_atomicnum)):
        for j in range(len(ALL_atomicnum)):
            ALL_1_atomicnum[i,j]=ALL_atomicnum[i]*ALL_atomicnum[j]
    #------------------------------#
    #------------------------------#
    #   Compute DISTANCE MATRIX    #
    #------------------------------#
    import scipy.spatial

    ALL_1 = scipy.spatial.distance_matrix(ALL_1, ALL_1)

    maximum_distance=np.amax(ALL_1)
    distance_matrix_uper=ALL_1[np.triu_indices(len(ALL_1),k=1)]
    atomicnum_uper=ALL_1_atomicnum[np.triu_indices(len(ALL_1_atomicnum),k=1)]
    occurrence, distances = np.histogram(distance_matrix_uper, bins = np.arange(0.0,maximum_distance,binn))
    occurrence_weighted, distances = np.histogram(distance_matrix_uper, bins = np.arange(0.0,maximum_distance,binn), weights=atomicnum_uper)
    # what's the highest distance in the array?
    print("Max distance in the 1 distance matrix: " + str(np.amax(ALL_1)))
    
    lines_occ=[]
    for i in range(len(occurrence)):
        lines_='{0: <50}'.format(distances[i])+'{0: <50}'.format(occurrence[i])
        lines_occ.append(lines_)
    distfilename='all_occur_'+str(snap)+'_'+str(binn)+'.txt'
    with open(distfilename, 'w') as f:
            f.write("#distance    occurrence\n")
            f.writelines("%s\n" % l for l in lines_occ)
    lines_occ=[]
    for i in range(len(occurrence_weighted)):
        lines_='{0: <50}'.format(distances[i])+'{0: <50}'.format(occurrence_weighted[i])
        lines_occ.append(lines_)
    distfilename='all_occur_weighted_'+str(snap)+'_'+str(binn)+'.txt'
    with open(distfilename, 'w') as f:
            f.write("#distance    occurrence\n")
            f.writelines("%s\n" % l for l in lines_occ)
    ####################################
    #  Computing the Fourier Transform #
    ####################################
    q_points = len(distances)
    #occurrence_1_ft = np.zeros(q_points, dtype=np.complex64)
    distance_ft = np.zeros(q_points, dtype=np.complex64)
    distance_ft_weighted = np.zeros(q_points, dtype=np.complex64)

    # Define "q_points" equally spaced q points in the interval 2*pi/d_max - 2*pi/d_min (d_max = 256 ang, d_min = 1 ang)
    qrange = np.arange(2*np.pi/q_points,2*np.pi/1,(2*np.pi/1-2*np.pi/q_points)/q_points)
    for q in range(0,len(qrange)):
        for ind_occ, occ  in enumerate(occurrence):
            if (distances[ind_occ]!=0.0 and qrange[q]>0.1):
                distance_ft[q]=distance_ft[q]+occ*np.sin(qrange[q]*distances[ind_occ])/qrange[q]/distances[ind_occ]
                distance_ft_weighted[q]=distance_ft_weighted[q]+occurrence_weighted[ind_occ]*np.sin(qrange[q]*distances[ind_occ])/qrange[q]/distances[ind_occ]
    
    lines_scat=[]
    for i in range(len(distance_ft)):
        lines_='{0: <50}'.format(qrange[i])+'{0: <50}'.format(len(occurrence)+2*np.real(distance_ft[i]))
        lines_scat.append(lines_)
    distfilename='all_scat_'+str(snap)+'_'+str(binn)+'.txt'
    with open(distfilename, 'w') as f:
            f.write("#qrange    normalized scattering\n")
            f.writelines("%s\n" % l for l in lines_scat)

    lines_scat=[]
    for i in range(len(distance_ft)):
        lines_='{0: <50}'.format(qrange[i])+'{0: <50}'.format(len(occurrence_weighted)+2*np.real(distance_ft_weighted[i]))
        lines_scat.append(lines_)
    distfilename='all_scat_weighted_'+str(snap)+'_'+str(binn)+'.txt'
    with open(distfilename, 'w') as f:
            f.write("#qrange    normalized scattering\n")
            f.writelines("%s\n" % l for l in lines_scat)



