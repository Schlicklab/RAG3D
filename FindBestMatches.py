import os, sys
from operator import itemgetter
from collections import OrderedDict

from RAG_v2 import Structure

#First C++ RAGTOP code must be run (It also removes pseudoknots.) 
#Loops, Vertices, VertexTypes, Graph, BPSEQ files are needed in the input directory.
#Requires subgraphs.py script which is the modified version of RNA matrix code. This version also creates subgraph info. It needs %d-vertex-modified files (eigenvalue info)
#PdbToNb.txt file is required for pdb number/pdb id conversion

#Under the output directory 3 subfolders are created to store the results: Results, Subgraphs, and Junctions
#This script also creates 2 subfolders: PML and MATCHES to store pymol script files & results of alignments, and the best matching results sorted wrt RMDSD values
#Best matches are searched against the RAG-3D database. The path to the database should be specified.

#S.J. 05/09/2018 - clean up the file     

#if len(sys.argv) != 4:
#   	print "Require 4 arguments: %s <QueryName> <Input Directory> <Output Directory>" %(sys.argv[0])
#   	sys.exit(1)

QueryName = None
InputDir = None
OutputDir = None
DatabaseDir = None
GraphName = None
matchVertexType = False # S.J. 05/10/2018 - to match vertex types only if flag specified

# S.J. 05/09/2018 - parse command line arguments
for i in range(1,len(sys.argv)):
	if(sys.argv[i] == "-query"):
		QueryName = sys.argv[i+1]
	elif(sys.argv[i] == "-input"):
		InputDir = sys.argv[i+1]
	elif(sys.argv[i] == "-output"):
		OutputDir = sys.argv[i+1]
	elif(sys.argv[i] == "-database"):
		DatabaseDir = sys.argv[i+1]
	elif(sys.argv[i] == "-graph"):
		GraphName = sys.argv[i+1]
   	elif(sys.argv[i] == "-matchVertexType"):
        	matchVertexType = True

if(QueryName == None):
	print "Query not specified"
	sys.exit(1)
if(InputDir == None):
	print "Input directory not specified"
	sys.exit(1)
if(OutputDir == None):
	print "Output directory not specified"
	sys.exit(1)
if(DatabaseDir == None):
	print "Database directory not specified"
	sys.exit(1)
if(GraphName == None):
	GraphName = ""

print "Name: %s, Graph file: %s_%sGraph,"%(QueryName,QueryName,GraphName),

#QueryName = sys.argv[1] #structure name without the extension
#InputDir = sys.argv[2]
#OutputDir = sys.argv[3]
#DatabaseDir = sys.argv[4]

try:
	os.makedirs(OutputDir)
except OSError:
	if os.path.exists(OutputDir):
		pass
	else:
		raise

output_path="%s"%OutputDir
PML_path="%sPML/"%OutputDir
MATCH_path="%sMATCHES/"%OutputDir


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

f1="%s%s.bpseq"%(InputDir,QueryName)
f2="%sLoops%s.txt"%(InputDir,QueryName)
f3="%s%s_%sGraph.pdb"%(InputDir,QueryName,GraphName)
f4="%sVertices%s.txt"%(InputDir,QueryName)
f5="File not needed"
f6="%sPdbToNb.txt"%(DatabaseDir)
if(matchVertexType): # S.J. 05/10/2018
    f7="%sVertexTypes%s.txt"%(InputDir,QueryName) # S.J. 05/07/2017 for reading in vertex types
else:
    f7="File not needed"

dbpath="%sResults/"%(DatabaseDir) #for search

fnc_details="%sdetails.txt"%(DatabaseDir) #this is for web tool. finds the fnc, experimental method and resolution data

R=Structure(QueryName,f1,f2,f3,f4,f5,f6,f7,output_path) # S.J. added f7
R.create_dirs()
R.readfiles()
#R.print_subgraphs()
R.calc_subgraphs()
R.find_matches(dbpath,fnc_details,matchVertexType) # S.J. 05/10/2018 adding the last argument
