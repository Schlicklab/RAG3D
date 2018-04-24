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

#S.J. 02/04/2016 - changing the command line arguments to take three arguments
if len(sys.argv) != 4:
   	print "Require 3 arguments: %s <QueryName> <Input Directory> <Output Directory>" %(sys.argv[0])
   	sys.exit(1)
    
QueryName = sys.argv[1] #structure name without the extension
#RandomNumber = sys.argv[2] #name output path
InputDir = sys.argv[2]
RandomNumber = sys.argv[3]

try:
	os.makedirs(RandomNumber)
except OSError:
	if os.path.exists(RandomNumber):
		pass
	else:
		raise

#output_path="./%s/"%RandomNumber
#PML_path="./%s/PML/"%RandomNumber
#MATCH_path="./%s/MATCHES/"%RandomNumber
output_path="%s"%RandomNumber
PML_path="%sPML/"%RandomNumber
MATCH_path="%sMATCHES/"%RandomNumber


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


#f1="./BPSEQ/%s.bpseq"%QueryName
#f2="./Loops-below-11v-numbers/Loops%s.txt"%QueryName
#f3="./GRAPH/%s_PDBGraph.pdb"%QueryName
#f4="./Vertices-below-11v-numbers/Vertices%s.txt"%QueryName
#f5="./PDB/%s.pdb"%QueryName
#f6="./PdbToNb.txt"

f1="%s%s.bpseq"%(InputDir,QueryName)
f2="%sLoops%s.txt"%(InputDir,QueryName)
f3="%s%s_PDBGraph.pdb"%(InputDir,QueryName)
f4="%sVertices%s.txt"%(InputDir,QueryName)
f5="%s%s.pdb"%(InputDir,QueryName)
f6="./PdbToNb.txt"


#dbpath="/ehome/cs4367/FullDataset/RAG-3D-26June/Results/" #for search
dbpath="/Users/sj78/sourcecodes/RAG3D_Database/Results/" #for search

fnc_details="./details.txt" #this is for web tool. finds the fnc, experimental method and resolution data

R=Structure(QueryName,f1,f2,f3,f4,f5,f6,output_path)
R.create_dirs()
R.readfiles()
#R.print_subgraphs()
R.calc_subgraphs()
R.find_matches(dbpath,fnc_details)
