#!/usr/local/bin/python
# Swati Jain (S.J.) - script for enumerating dual graph using the rules set in various papers
# 07/06/2018 - takes three integers as command line arguments, 1: number of vertices for which graphs have to be enumerated, 2: number of vertices of subgraph 1, 3: number of vertices of subgraph 2. All graphs with sub 1 vertices and all graphs with sub 2 vertices will be combined to generate new graphs. Assumption: sub 1 vertex >= sub 2 vertex, therefore sub 2 graphs will be added to sub 1 graph

import sys
from copy import deepcopy
import numpy.linalg as LA
from decimal import *
import os
import time
from numpy import *
from ClassesFunctions import *

def readSubTreeGraphs(vertex,Graphs):

    prevAdjFile = "V%dAdjDG"%(vertex) # output files
    prevEigenFile = "%dEigen"%(vertex)

    loadEigenvalues(Graphs,vertex,prevEigenFile)
    loadAdjMatrices(Graphs,vertex,prevAdjFile) #reading the adjacency matrices

# function to determine if the laplacian of the dual graph has already been enumerated or not - 05/16/2018
def determineUnique(adjMatrix):
    
    global total_new_graphs
    # calculate the laplacian and eigen values
    decimalArray = calcEigenValues(adjMatrix)
    eigen_tuple=tuple(decimalArray) # 09/14/2018 - implementation of the search using a python dictionary

    # 09/14/2018 - implementation of the search using a python dictionary
    Graphs=UniqueGraphs.get(eigen_tuple,[])
    flag=False
    for g in Graphs:
	    flag=checkIsomorphism(g.adjMatrix,adjMatrix)
            if flag: # if this is already present, just return
                #return
                break

    # code will come here if the current graph is not present, so needs to be added
    if not flag:
        total_new_graphs+=1
        new_g=TreeGraph(len(adjMatrix),adjMatrix,total_new_graphs,decimalArray)
        new_g.printEigen(eigenFile)
        printMat(new_g.adjMatrix,adjMatFile)
        Graphs.append(new_g)
        UniqueGraphs[eigen_tuple]=Graphs

# 09/17/2018 - generate graphs for each pair of g1 and g2, so that only edge combinations that are compatible are generated
def enumerateAll():

    basicAdjMat = [[0] * (numVertices) for a in range(numVertices)] # the starting matrix with every element 0

    for g in Sub1_Graphs: # for each sub graph to which one vertex has to be added
        
        basicAdjMat[0:numVertices] = [[0] * (numVertices) for a in range(numVertices)]
        for i in range(0,numVertices-1): # copying the subgraph adjMatrix
            for j in range(0,numVertices-1):
            	basicAdjMat[i][j] = g.adjMatrix[i][j]
        
	for vert in range(0,numVertices-1):

	    curAdjMatrix = deepcopy(basicAdjMat)
	    curAdjMatrix[vert][numVertices-1] = 1
	    curAdjMatrix[numVertices-1][vert] = 1

            determineUnique(curAdjMatrix)


# The main function
numVertices = int(sys.argv[1])

UniqueGraphs = dict() # 09/14/2018 - implementation of the search using a python dictionary
total_new_graphs=0

Sub1_Graphs = [] 

adjMatFile = "V%dAdjDG"%(numVertices) # output files for newly enumerated graphs
eigenFile = "%dEigen"%(numVertices)

readSubTreeGraphs(numVertices-1,Sub1_Graphs) # read the eigen values and adjacency matrices for dual graphs for Sub1

enumerateAll()
