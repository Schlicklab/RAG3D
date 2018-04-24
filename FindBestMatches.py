import os, sys
from operator import itemgetter
from collections import OrderedDict

from RAG_v2 import Structure

#First C++ code must be runned. (It also removes pseudoknots.) 
#Loops, Vertices, PDBGraph, BPSEQ files are needed.
#Requires subgraphs.py script which is the modified version of RNA matrix code. This version also creates subgraph info. It needs %d-vertex-modified files (eigenvalue info)
#PdbToNb.txt file is required for pdb number/pdb id conversion

#The output files are written under directory named <RandomNumber>.
#Under this directory 3 subfolders are created to store the results: Results, Subgraphs, and Junctions
#This script also creates 2 subfolders: PML and MATCHES to store pymol script files & results of alignments, and the best matching results sorted wrt RMDSD values
#Best matches are searched against the RAG-3D database. The path to the database should be specified.


if len(sys.argv) != 3:
   	print "Require 2 arguments: %s <QueryName> <RandomNumber>" %(sys.argv[0])
   	sys.exit(1)
    
QueryName = sys.argv[1] #structure name without the extension
RandomNumber = sys.argv[2] #name output path

try:
	os.makedirs(RandomNumber)
except OSError:
	if os.path.exists(RandomNumber):
		pass
	else:
		raise

output_path="./%s/"%RandomNumber
PML_path="./%s/PML/"%RandomNumber
MATCH_path="./%s/MATCHES/"%RandomNumber


try:
	os.makedirs(PML_path)
except OSError:
	if os.path.exists(PML_path):
		pass
	else:
		raise

try:
	os.makedirs(MATCH_path)
except OSError:
	if os.path.exists(MATCH_path):
		pass
	else:
		raise


#f1="../FullDataset/BPSEQ/%s.bpseq"%QueryName
#f2="../FullDataset/Loops-below-11v-numbers/Loops%s.txt"%QueryName
#f3="../FullDataset/GRAPH/%s_PDBGraph.pdb"%QueryName
#f4="../FullDataset/Vertices-below-11v-numbers/Vertices%s.txt"%QueryName
#f5="../FullDataset/PDB/%s.pdb"%QueryName
f6="./PdbToNb.txt"
f1="../3IYR/%s.bpseq"%QueryName
f2="../3IYR/Loops%s.txt"%QueryName
f3="../3IYR/%s_PDBGraph.pdb"%QueryName
f4="../3IYR/Vertices%s.txt"%QueryName
f5="../3IYR/%s.pdb"%QueryName
	
#dbpath="/ehome/cs4367/FullDataset/RAG-3D-26June/Results/" #for search
dbpath="/ehome/cs4367/FullDataset/RAG-3D/Results/" #for search

fnc_details="./details.txt" #this is for web tool. finds the fnc, experimental method and resolution data

R=Structure(QueryName,f1,f2,f3,f4,f5,f6,output_path)
R.create_dirs()
R.readfiles()
#R.print_subgraphs()
R.calc_subgraphs()
R.find_matches(dbpath,fnc_details)
