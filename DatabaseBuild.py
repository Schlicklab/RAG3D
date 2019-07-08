import os, sys
from operator import itemgetter
from collections import OrderedDict

from RAG_v2 import Structure

#It reads the information previosly created with the C++ code (Loops,Vertices,Graph)
#Requires subgraphs_naoto.py script which is the modified version of RNA matrix code. This version also creates subgraph info. It needs %d-vertex-modified files (eigenvalue info)
#PdbToNb.txt file is required for pdb number/pdb id conversion or just a list of pdb id/pdb id if not converted to numbers
#Under output_path 3 directories are created: Results, Subgraphs, Junctions
# S.J. - 05/07/2019 - please update the paths when using this file
# sys.argv[1] is a list of PDB files to be processed (without extension)

output_path="/ts_home/sj78/labwork/RAG-3D-Database/"
#output_path="./RAG-3D-Test/"

try:
	os.makedirs(output_path)
except OSError:
	if os.path.exists(output_path):
		pass
	else:
		raise
list=open(sys.argv[1])
files = list.readlines()
list.close()
for pdb in files:
	pdb = pdb.strip()
	print pdb
    
    	f1="/Users/sj78/Documents/labwork/Amiel/FinalBPSEQs_LC_noPK/%s.bpseq"%pdb
    	f2="/Users/sj78/Documents/labwork/Amiel/FinalBPSEQs_LC_noPK/LoopsVertices/Loops/Loops%s.txt"%pdb
    	f3="/Users/sj78/Documents/labwork/Amiel/FinalBPSEQs_LC_noPK/LoopsVertices/Graphs/%s_PDBGraph.pdb"%pdb
    	f4="/Users/sj78/Documents/labwork/Amiel/FinalBPSEQs_LC_noPK/LoopsVertices/Vertices/Vertices%s.txt"%pdb
    	f5="/ts_home/sj78/labwork/AllRNADataset_August2018/RNAPDBs/%s.pdb"%pdb
    	f6="/ts_home/sj78/labwork/RAG-3D-Database/PdbToNb.txt"
	f7="/Users/sj78/Documents/labwork/Amiel/FinalBPSEQs_LC_noPK/LoopsVertices/VertexTypes/VertexTypes%s.txt"%pdb

	#S.J.R=Structure(pdb,f1,f2,f3,f4,f5,f6,output_path)
	R=Structure(pdb,f1,f2,f3,f4,f5,f6,f7,output_path)
	R.create_dirs()
	if R.readfiles()==False:
		continue
	else:
		R.calc_subgraphs()
