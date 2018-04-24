import os, sys
from operator import itemgetter
from collections import OrderedDict

from RAG_v2 import Structure

#It reads the information previosly created with the C++ code (Loops,Vertices,Graph)
#Requires subgraphs.py script which is the modified version of RNA matrix code. This version also creates subgraph info. It needs %d-vertex-modified files (eigenvalue info)
#PdbToNb.txt file is required for pdb number/pdb id conversion
#Under output_path 3 directories are created: Results, Subgraphs, Junctions

output_path="../FullDataset/RAG-3D/"
try:
	os.makedirs(output_path)
except OSError:
	if os.path.exists(output_path):
		pass
	else:
		raise

exception=[366,2299,1975,2214,2215,2216,2217] #not working pdb numbers

for i in range(1,2327): #for each pdb number
#for i in range(936,937): #for each pdb number
	if i in exception:
		continue

	pdb=str(i)
	f1="../FullDataset/BPSEQ/%s.bpseq"%pdb
	f2="../FullDataset/Coutput/Loops%s.txt"%pdb
	f3="../FullDataset/GRAPH/%s_PDBGraph.pdb"%pdb
	f4="../FullDataset/Coutput/Vertices%s.txt"%pdb
	f5="../FullDataset/PDB/%s.pdb"%pdb
	f6="./PdbToNb.txt"

#large=os.listdir("../data-03-29-14/bpseq-above-10v/")
#for i in large: #for each pdb number
#
#	pdb=i.split('.')[0]
#	print pdb
#	f1="../data-03-29-14/bpseq-above-10v/%s.bpseq"%pdb
#	f2="../data-03-29-14/Loops-above-10v/Loops%s.txt"%pdb
#	f3="../data-03-29-14/PDBGraph-above-10v/%s_PDBGraph.pdb"%pdb
#	f4="../data-03-29-14/Vertices-above-10v/Vertices%s.txt"%pdb
#	f5="../data-03-29-14/pdb-above-10v/%s.pdb"%pdb
#	f6="./PdbToNb.txt"

	R=Structure(pdb,f1,f2,f3,f4,f5,f6,output_path)
	R=Structure(pdb,f1,f2,f3,f4,f5,f6,output_path)
	R.create_dirs()
	if R.readfiles()==False:
		continue
	else:
		R.calc_subgraphs()
