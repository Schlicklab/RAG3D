import os, sys
from operator import itemgetter
from collections import OrderedDict
import numpy
import time

class Structure:
	def __init__(self,name,bpseqfile,loopfile, graphfile, vertfile, pdbfile, pdbnbfile, vertTypefile, outpath):
		self.name=name 
		self.bpseqfile=bpseqfile
		self.loopfile=loopfile
		self.graphfile=graphfile
		self.vertfile=vertfile
		self.pdbfile=pdbfile
		self.pdbnbfile=pdbnbfile
        	self.vertTypefile=vertTypefile # S.J. 05/07/2017
		self.outpath=outpath
		
		#Under outpath create 3 directories
		self.Results_dir=self.outpath + "Results/" #includes all graph IDs as directories. GRAPH-FRGT and AllAtom are subdirectories.
		self.Subgraphs_dir=self.outpath + "Subgraphs/" #includes results of subgraphs.py i.e., subgraph IDs and corresponding vertices info.
		self.Junctions_dir=self.outpath + "Junctions/" #includes 3way and 4way junction info.

		self.Pdb_Nb={} #pdb id to pdb number conversion 
		self.loopinfo=[] 
		self.graphinfo=[]
		self.pdbinfo=[]
		self.vertinfo=[]
		self.vertTypeinfo=[] # S.J. 05/09/2017
		self.ID={} #subgraph id and corresponding vertices. 2_1-1:[0,1] or 3_1-2:[0,1,2]
		self.sorted_allinfo=[] #loop & graph files info.
		self.Vert={} #ordering number and ibp from loops file & atom no from graph file. 0:[1,3,18,28,1]
		
		self.GraphID= 0

	def create_dirs(self):
		##Results contains subfolders named after graph ids (2_1, 3_1, etc.)
		##Under each subfolder there are GRAPH-FRGT and AllAtom folders containing corresponding pdb files
		try:
			os.makedirs(self.Results_dir)
		except OSError:
			if os.path.exists(self.Results_dir):
				pass
			else:
				raise

		##Subgraphs contains output files of subgraph.py script
		try:
			os.makedirs(self.Subgraphs_dir)
		except OSError:
			if os.path.exists(self.Subgraphs_dir):
				pass
			else:
				raise

		##Junctions include 3-way and 4-way junction info
		try:
			os.makedirs(self.Junctions_dir)
		except OSError:
			if os.path.exists(self.Junctions_dir):
				pass
			else:
				raise


	def readfiles(self):
		
		#subgraph.py is the modified version of the RNA matrix code
		sfile=self.Subgraphs_dir + self.name + "-subgraphs"
		os.system('/usr/bin/python2.7 subgraphs_naoto.py  %s | grep "Subgraph" | sort -nrk 3,3 > %s'%(self.bpseqfile,sfile))

		try:
			file1=open(sfile,"r")
		except IOError:
			print >> sys.stderr, "Subgraphs file could not be opened"
			return False
			#sys.exit()
		sublines=file1.readlines()
		if len(sublines)==0:
			print >> sys.stderr, "Subgraphs file is empty"
			return False
			#sys.exit()
		file1.close()

		try:
			file2=open(self.loopfile,"r")
		except IOError:
			print >> sys.stderr, "Loops file could not be opened"
			return False
			#sys.exit()
		looplines=file2.readlines()
		file2.close()

		try:
			file3=open(self.graphfile,"r")
		except IOError:
			print >> sys.stderr, "PDB Graph file could not be opened"
			return False
			#sys.exit()
		graphlines=file3.readlines()
		file3.close()

		try:
			file4=open(self.vertfile,"r")
		except IOError:
			print >> sys.stderr, "Vertices file could not be opened"
			#sys.exit()
			return False
		vertlines=file4.readlines()
		file4.close()

		#try: S.J. commented out to not read the pdb file
		#	file5=open(self.pdbfile,"r")
		#except IOError:
		#	print >> sys.stderr, "Pdb file could not be opened"
		#	#sys.exit()
		#	return False
		#pdblines=file5.readlines()
		#file5.close()

		try:
			file6=open(self.pdbnbfile,"r")
		except IOError:
			print >> sys.stderr, "PdbToNb file could not be opened"
			return False
			#sys.exit()
		nblines=file6.readlines()
		file6.close()

        	# S.J. 05/07/2017 - added to read the vertType file
        	try:
            		file7=open(self.vertTypefile,"r")
        	except IOError:
            		print >> sys.stderr, "Verttype file could not be opened"
            		return False
            		#sys.exit()
        	vertTypelines=file7.readlines()
        	file7.close()

		try:
			infofile=open(self.outpath + "Query_info.txt","w") #for web page
		except IOError:
			print >> sys.stderr, "Query info file could not be opened"

		#Write into arrays or dictionaries
		for line in nblines:
			record=line.split()
			self.Pdb_Nb.update({record[1]:record[0]}) #PdbNumber: PdbID_ChainID
		#print self.Pdb_Nb

		AllNgraph=[]
		for i in range(len(graphlines)):
			self.graphinfo.append(graphlines[i])
			if 'ATOM' in graphlines[i]:
				columns= graphlines[i].split()
				if columns[2].strip()=="N":
					AllNgraph.append(graphlines[i])
		Ngraph=[]
		all_loops=[]
		for i in range(len(looplines)):
			if "(loop" in looplines[i]:
				all_loops.append(looplines[i])
			elif "3-way" in looplines[i]:
				columns=looplines[i].split()
				num_of_3way=int(columns[3]) #number of 3-way junctions
				j1=open(self.Junctions_dir+ self.name +"-3way.txt","w")
				for x in range(num_of_3way+1):
					j1.write("%s"%looplines[i+x])
				j1.close()
			elif "4-way" in looplines[i]:
				columns=looplines[i].split()
				num_of_4way=int(columns[3]) #number of 4-way junctions
				j2=open(self.Junctions_dir + self.name + "-4way.txt","w")
				for x in range(num_of_4way+1):
					j2.write("%s"%looplines[i+x])
				j2.close()
		#print self.loopinfo

		ignored_loops=[]
		num_of_internal=0
		num_of_hairpin=0
		for i in range(len(all_loops)):
			if "Internalloop" in all_loops[i]:
				num_of_internal+=1
				columns=all_loops[i].split()
				if int(columns[4])-int(columns[2])<=2 and int(columns[8])-int(columns[6])<=2:
					temp=[]
					temp.append(int(columns[2]))
					temp.append(int(columns[4]))
					ignored_loops.append(temp)
					temp=[]
					temp.append(int(columns[6]))
					temp.append(int(columns[8]))
					ignored_loops.append(temp)
					continue
				else:
					self.loopinfo.append(all_loops[i])
					Ngraph.append(AllNgraph[i])
			else:
				self.loopinfo.append(all_loops[i])
				Ngraph.append(AllNgraph[i])
				if "Hairpin" in all_loops[i]:
					num_of_hairpin+=1
		for item in self.loopinfo:
			infofile.write(item)
		infofile.write("\nTotal Hairpins: %d\n"%num_of_hairpin)
		infofile.write("Total Internal loops: %d\n"%num_of_internal)
		try:
			num_of_3way
			infofile.write("Total 3-way junctions: %d\n"%num_of_3way)
		except NameError:
			junc_exists=False
		try:
			num_of_4way
			infofile.write("Total 4-way junctions: %d\n"%num_of_4way)
		except NameError:
			junc_exists=False

		#S.J. commented out to not use the pdb file
		#for i in range(len(pdblines)):
		#	if 'ATOM' in pdblines[i]:
		#		self.pdbinfo.append(pdblines[i])
		#print self.pdbinfo
		for i in range(len(vertlines)):
           		if vertlines[i].startswith("Vertices:"):
				check=0
				for out in ignored_loops:
					if int(vertlines[i].split()[3]) in out:
						check=1
				if check==0:
        				self.vertinfo.append(map(int,[vertlines[i].split()[1],vertlines[i].split()[3]]))
		#print "vertices info", self.vertinfo
	
    
		for i in range(len(vertTypelines)):
			x=vertTypelines[i].split(":")[0]
			y=x.split()[0]
			z=int(y)
			st=vertTypelines[i].split(":")[1]
			st1=st.split("\n")[0]
			self.vertTypeinfo.append([z,st1])
		#print "vert type info", self.vertTypeinfo[0]
    
		IDlist=[]	
		for i in range(len(sublines)):
			columns=sublines[i].split()
			sub_ID=columns[2]
			vertices=map(int,columns[3:])
			IDlist.append([sub_ID,vertices])

		self.GraphID=IDlist[0][0]
		print "Name: %s, Graph ID: %s"%(self.name,self.GraphID)
		infofile.write("\nTotal number of vertices: %s\n"%self.GraphID.split('_')[0])
		infofile.write("\nQuery has the topology ID: %s\n"%self.GraphID)


		#add numbers at the end of the subgraph IDs. For example: 2_1-1, 2_1-2,...
		for i in range(len(IDlist)):
			k=IDlist[i][0]+'-'+str(i)
			v=IDlist[i][1]
			self.ID.update({k:v})		

		infofile.write("\nQuery can be divided into %d subgraphs\n"%(len(self.ID.items())-1))
		allinfo=[] #save base indices and atom numbers as a list of lists
		for i in range(len(self.loopinfo)):
			columns=self.loopinfo[i].split()
			temp=[]
			atomno=int(Ngraph[i].split()[1])
			for j in range(0,len(columns)-2,2):
				temp.extend([int(columns[j+2])])
			temp.extend([atomno])
			allinfo.append(temp)


		self.sorted_allinfo=sorted(allinfo, key=lambda x:int(x[0])) #sort wrt first column

		for i in range(len(self.sorted_allinfo)):
			self.Vert.update({i:self.sorted_allinfo[i]})
		#print self.sorted_allinfo
		#print self.Vert

		return True
            
	def print_subgraphs(self):
		#### Print Subgraph id and vertices
		print "Subgraphs:"
		for keys,values in self.ID.items():
			print keys, values
		print "\n"

	def calc_subgraphs(self):
		if len(self.loopinfo)!=int(self.GraphID.split('_')[0]):
			print "Number of vertices does not match!"
			#sys.exit() 
		else:
			for key in self.ID:
				loops_in=[] #loops involved in the subgraph
				loops_out=[] #loops not involved in the subgraph
				graph=[]
				forAAextract=[]
				for v in self.Vert.keys():
					pair=zip(self.Vert[v][::2],self.Vert[v][1::2]) #read each loop interval as a pair
					if v in self.ID[key]:
						loops_in.extend(set(pair)) #remove duplicates by using set
						graph.append(self.Vert[v][-1]) #append N atom numbers to graph file
					else:
						loops_out.extend(set(pair))
				sorted_loop_in= sorted(loops_in, key=lambda t: int(t[0]))
				sorted_loop_out= sorted(loops_out, key=lambda t: int(t[0]))
				if len(sorted_loop_out)<1:
					mini=sorted_loop_in[0][0]
					maxi=sorted_loop_in[-1][1]
					forAAextract.append((mini,maxi+1))
				else:
				
					offset1=0 #for iterating loop_in
					offset2=0 #for iterating loop_out
					for i in range(len(sorted_loop_out)):
						if sorted_loop_out[i][0]<sorted_loop_in[0][0]: #if the smallest bp is in loop_out
							offset2=i+1
						else:
							break
					flag=True
					while(flag):
						
						mini=sorted_loop_in[offset1][0]
						maxi=sorted_loop_in[offset1][1]
						for i in range(offset1,len(sorted_loop_in)):
							if offset2==len(sorted_loop_out) or sorted_loop_in[i][1]<sorted_loop_out[offset2][0]:
								maxi=sorted_loop_in[i][1]
							else:
								
								forAAextract.append((mini,maxi+1))
								offset1=i
								for j in range(0,len(sorted_loop_out)):
									offset2+=1
									if offset2==len(sorted_loop_out):
										break
									if sorted_loop_out[offset2][0]>sorted_loop_in[offset1][0]:
										break
								break
							if i==len(sorted_loop_in)-1:
								maxi=sorted_loop_in[i][1]
								forAAextract.append((mini,maxi+1))
								flag=False
								break
	
				#append C atom numbers to graph file
				for v in self.vertinfo:
					for i in range(len(forAAextract)):
						if v[1] in range(*forAAextract[i]):
							if v[0] not in graph:
								graph.append(v[0])
				forGraph=sorted(graph)
				path1=self.Results_dir + key.split('-')[0] + "/AllAtom/"
				name1="%s-%s-AA-frgt.pdb"%(self.name, key)
				self.write_AA(path1,name1,forAAextract)

				path2=self.Results_dir + key.split('-')[0] + "/GRAPH-FRGT/"
				name2="%s-%s-3D-frgt.pdb"%(self.name, key)
				self.write_graph(path2,name2,forGraph)
			

	def write_AA(self,AA_path,filename,values):
		
		#AA_path=self.Results_dir + self.GraphID + "/AllAtom/"
		try:
			os.makedirs(AA_path)
		except OSError:
			if os.path.exists(AA_path):
				pass
			else:
				raise
		try:
			frgt=open(AA_path + filename,"w")
		except IOError:
			print >> sys.stderr, "AA-frgt file could not be opened"
			sys.exit()

		
		#diff=int(self.pdbinfo[0][22:26])-1 #read the number of first atom. for in case if it does not start with 1! the difference btw atom no in pdb and bpseq
		for line in self.pdbinfo:
			#atom=int(line[22:26])
			atom=1
			for i in range(len(values)):
				#myrange=(values[i][0]+diff,values[i][1]+diff) #bpseq starts with 1. if pdb atom no starts with >1, add the difference!
				myrange=(values[i][0],values[i][1])
				if atom in range(*myrange): #star is used to unpack the tuple (x,y)
					frgt.write("%s"%line)
			atom+=1
		frgt.write("END\n")
		frgt.close()

	def write_graph(self,Graph_path,filename,values):
		#Graph_path=self.Results_dir + self.GraphID + "/GRAPH-FRGT/"
		try:
			os.makedirs(Graph_path)
		except OSError:
			if os.path.exists(Graph_path):
				pass
			else:
				raise
		try:
			frgt3d=open(Graph_path + filename,"w")
		except IOError:
			print >> sys.stderr, "3d-frgt file could not be opened"
			sys.exit()
		mylines=[]
		atoms=[]
		connections=[]
		##Read graph pdb file
		for line in self.graphinfo:
			columns= line.split()
			if columns[0]=='ATOM':
				atoms.append(int(columns[1])) #save all atom numbers
			elif columns[0]=='CONECT':
				temp=[]
				temp.append(int(columns[1]))
				temp.append(int(columns[2]))
				connections.append(temp)	#save all connections as pairs
	
		out_val=[] #vertices those not included in the considered subgraph
		for a in atoms:
			if a not in values:
				out_val.append(a)
		for v in out_val: #rearrange connections for the broken vertices
			size=len(connections)
			to_remove=[]
			old_pair=[]
			for i in range(size):
				if v==int(connections[i][0]):
					if connections[i][1] not in old_pair: #and connections[i][1] not in out_val:
						old_pair.append(connections[i][1])
						to_remove.append(connections[i])
				elif v==int(connections[i][1]):
					if connections[i][0] not in old_pair: #and connections[i][0] not in out_val:
						old_pair.append(connections[i][0])
						to_remove.append(connections[i])
			if len(old_pair)==2:
				connections.append(old_pair) #connect the broken vertices

			elif len(old_pair)>2:#for the nodes with more than 2 connections
				for x in range(1,len(old_pair)-1):
					connections.append([old_pair[0],old_pair[x]]) #connect the broken vertices
			for i in range(len(to_remove)):
				connections.remove(to_remove[i]) #remove connections
		#print connections
		##Write new graph PDB file
		for line in self.graphinfo:
			columns= line.split()
			if columns[0]=='ATOM':
				if int(columns[1]) in values:
					frgt3d.write("%s"%line)
			elif columns[0]=='TER':
				frgt3d.write("%s"%line)
		
		for i in range(len(connections)):
			frgt3d.write("%-6s%5d%5d\n"%("CONECT",int(connections[i][0]),int(connections[i][1])))
				
		frgt3d.close()
	
	# S.J. 05/09/2017 function added
	def isSimilar(self,key,dbpath,fragment):
		target_info=[]
		fragment_info=[]
		
		file=open("%s/Results/%s/GRAPH-FRGT/%s-%s-3D-frgt.pdb"%(self.outpath,key.split('-')[0],self.name,key),"r")
		targetlines=file.readlines()
		file.close()
		#getting vertex type for vertices in the target graph
		for i in range(len(targetlines)):
			if targetlines[i].split()[0] == "ATOM":
				if targetlines[i].split()[2] == "N":
					vertex=int(targetlines[i].split()[5])
					for j in range(len(self.vertTypeinfo)):
						if self.vertTypeinfo[j][0] == vertex:	
							target_info.append(self.vertTypeinfo[j][1])
		sorted_target=sorted(target_info)
		#print sorted_target
		
		#getting vertex type for vertices in the fragment graph
		fragmentpdb=fragment[3].split('-')[0] # number of the pdb id from which the fragment is
		file=open("%sVertexTypes/VertexTypes%s.txt"%(dbpath,fragmentpdb),"r")
		fragmentlines=file.readlines()
		file.close()
		fragment_full_info=[]
		for i in range(len(fragmentlines)): # reading the full pdb vertex type from which the fragment comes
			x=fragmentlines[i].split(":")[0]
                        y=x.split()[0]
                        z=int(y)
                        st=fragmentlines[i].split(":")[1]
                        st1=st.split("\n")[0]
                        fragment_full_info.append([z,st1])
		
		file=open("%s%s/GRAPH-FRGT/%s-3D-frgt.pdb"%(dbpath,key.split('-')[0],fragment[3]),"r")
		fragmentlines=file.readlines()
		file.close()
		for i in range(len(fragmentlines)):
			if fragmentlines[i].split()[0] == "ATOM":
				if fragmentlines[i].split()[2] == "N":
					vertex=int(fragmentlines[i].split()[5])
					for j in range(len(fragment_full_info)):
						if fragment_full_info[j][0] == vertex:
							fragment_info.append(fragment_full_info[j][1])
		sorted_fragment=sorted(fragment_info)
		#print sorted_fragment
		
		# check if the fragment and target have the same vertex types
		flag=1
		for i in range(len(sorted_target)):
			if sorted_target[i] == sorted_fragment[i]:
				flag=1
			else:
				flag=0
				break
		
		return flag 

	def find_matches(self,dbpath,inputfile):
		time1=time.time()
		fncfile=open(inputfile,"r")
		fnclines=fncfile.readlines()
		fncfile.close()
		ordered_ID=OrderedDict(sorted(self.ID.items(), key=lambda t: int(t[0].split('-')[1])))
		allmatches=open("%s/MATCHES/allmatches"%(self.outpath),"w")
		for key in ordered_ID:
			if os.path.isdir(dbpath + "%s/"%key.split('-')[0]):
				mypath=dbpath + "%s/GRAPH-FRGT/"%key.split('-')[0]
				aapath=dbpath + "%s/AllAtom/"%key.split('-')[0]
				fragments= os.listdir(mypath)	
				allresults=[]
				pml=open("%s/PML/%s.pml"%(self.outpath,key),"w")
				pml.write("load %s/%s-%s-3D-frgt.pdb\n"%(self.Results_dir+key.split('-')[0]+"/GRAPH-FRGT",self.name,key))
				for f in fragments:
					results=[]
					if self.name in self.Pdb_Nb:
						results.append(self.Pdb_Nb[self.name]) #pdb id of the query (if exists in our conversion table!)
					else:
						results.append(self.name)
					order=str(int(key.split('-')[1])+1)
					results.append(key.split('-')[0]+' ('+order+')') #graph id of the query
					results.append(self.Pdb_Nb[f.split('-')[0]]) #pdb id of the fragment
					results.append('-'.join(f.split('-')[:3])) #subgraph id of the fragment
					allresults.append(results)
					pml.write("load %s/%s\n"%(mypath,f))
					pml.write("align %s, %s-%s-3D-frgt\n"%(f[0:-4],self.name,key))
				pml.close()
				print "Finding matches for %s...\n"%key
				os.system("/Applications/MacPyMOL.app/Contents/MacOS/MacPyMOL -cq %s/PML/%s.pml > %s/PML/Alignment-%s.out"%(self.outpath,key,self.outpath,key))
				resultsfile=open("%s/PML/Alignment-%s.out"%(self.outpath,key),"r")
				mylines=resultsfile.readlines()
				resultsfile.close()
				j=0
				for line in mylines:
					columns=line.split()
					if columns[1]=="%s-%s-3D-frgt,"%(self.name,key): #this if statement added by S.J. - 01/20/2016 - as the rmsd between the same structure was not being calculated
						rmsd=0.000
						allresults[j].append(rmsd)
						j+=1
					if columns[0]=="Executive:":
						rmsd=float(columns[3])
						allresults[j].append(rmsd)
						j+=1
				sorted_results=sorted(allresults, key=lambda x:float(x[4])) #sort wrt rmsd
				#sorted_results=allresults
				#for i in range(len(sorted_results)):
					#print "Sorted: %s"%(sorted_results[i])

				best_results=[]#10 best hits
				count=0
				for i in range(len(sorted_results)-1):
					if count<10:
						best_results.append(sorted_results[i]) 
						count+=1
					if i==len(sorted_results)-1:
						break
				#tenmatches=open("%s/MATCHES/%s-10matches"%(self.outpath,key),"w")
				#for i in range(len(best_results)):
				order=0;
				for i in range(len(sorted_results)):
					similar=self.isSimilar(key,dbpath,sorted_results[i]) # S.J. 05/09/2017 for checking if the fragment 
					if similar==0:
						continue
					#print "Best Results Before=%s"%(best_results[i][2])
					#os.system("cp %s/%s-3D-frgt.pdb %s/MATCHES/"%(mypath,best_results[i][3],self.outpath))
					#os.system("cp %s/%s-AA-frgt.pdb %s/MATCHES/"%(aapath,best_results[i][3],self.outpath))
					for line in fnclines:
						columns=line.split('\t')
						if sorted_results[i][2]==columns[0]:
							os.system("cp %s/%s-3D-frgt.pdb %s/MATCHES/"%(mypath,sorted_results[i][3],self.outpath))
							os.system("cp %s/%s-AA-frgt.pdb %s/MATCHES/"%(aapath,sorted_results[i][3],self.outpath))
							frag=sorted_results[i][3]
							top=sorted_results[i][1]
							order+=1
							#order=i+1
							pdb=sorted_results[i][2]
							fnc=columns[1]
							method=columns[2]
							res=columns[3].strip()
							rmsd=sorted_results[i][4]
							#rmsd=0
							#print top, order, pdb, fnc, method, res, rmsd
							#tenmatches.write("%s\t%s\t%d\t%s\t%s\t%s\t%s\t%.3f\n"%(frag,top,order,pdb,fnc,method,res,rmsd))
							allmatches.write("%s\t%s\t%d\t%s\t%s\t%s\t%s\t%.3f\n"%(frag,top,order,pdb,fnc,method,res,rmsd))
							break
					if order==10:
						break
				#tenmatches.close()
				allmatches.write('\n')
			else:
				continue
		allmatches.close()
		time2=time.time()
		print "Time: %f seconds\n"%(float(time2-time1)/60)	
