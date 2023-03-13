# RAG3D

RAG-3D is a dataset of RNA tertiary (3D) structures and substructures designed to exploit graph representations of RNAs for the goal of searching for similar 3D structural fragments. The objects in RAG-3D consist of 3D structures translated into 3D graphs, cateloged based on the connectivity between their secondary structure elements. Each graph is additionally described in terms of its subgraph building blocks. The graph-substructuring approach for RNA structural search and classification is based on a hierarchical classification of RNA folds. These folds are represented as coarse-grained, tree-graph models, where RNA double-helices are represented as edges and loop domains (hairpins, internal loops, junctions and helices ends) are denoted as vertices. The RAG-3D search tool then compares a query RNA 3D structure to those in the database to obtain structurally similar structures and substructures. Since tree graphs cannot handle pseudoknots at present, pseudoknots are removed both from the query structure and from the structures in RAG-3D following the Elimination Gain (EG) method (Smit et al. 2008). The RAG-3D flowchart is shown in  ![Figure 4 in Zahran et al. NAR 2015](https://academic.oup.com/view-large/figure/81899349/gkv823fig4.jpeg). We implemented a user friendly web server to search our [RAG-3D database](http://www.biomath.nyu.edu/?q=RAG3D) for structure and substructure similarities. 


## From 2D graph to 3D graph representation
Each RNA 2D structure has a 2D tree graph representation based on the connectivity of its secondary motifs and is associated with a Laplacian matrix **L=D-A** which describes its connectivity. We expand RAD 2D into 3D graphs, as developed in (Laing et al. 2013, Kim et al. 2014). In RAG-3D tree graphs, the vertices are points in 3D space and the edges are line segments between vertices. We further incorporate junction feature details to define '3D' graphs as detailed previously (Kim et al. 2014, Kim et al. 2015). Verticces are added at the terminal base pairs of a helix to represent helices of different lengths. A central vertex is also added at the center of the junction domain to capture the junction's spatial properties. We also define additional edges to connect vertices at the end of helices, and to connect the center of the junction. 

We represent each helix by 2 vertices and 1 edge. The 3D coordinates of each vertex are determined in three steps: 

(i) find the midpoint M of C1′ atoms between the purine (Adenine and Guanine) and pyrimidine (Cytosine and Uracil) of the terminal base pairs of a helix;

(ii) consider the orthogonal projection from M on to the line connecting the C8 and C6 atoms of the purine and pyrimidine, respectively;

(iii) scale the vector projection by 4 Å. 

## Graph partitioning/subgraphs extraction

The initial graph partitioning is performed at the 2D level. The main condition guiding our 2D graph partitioning/subgraph identification protocol is Junction Intactness: a subgraph containing a junction must retain all its connected vertices. Our algorithm for subgraph extraction starts by determining from the diagonal matrix, D, whether the query structure contains junctions, i.e. vertices with more than two connections. The second step identifies all non-zero elements of each line of the adjacency matriix, A, to constitute an initial set of subgraphs. The third step is to identify all possible subgraphs, by merging the sets that share at least one vertex. Once we extract all subgraphs, we calculate a new Laplacian matrix for each subgraph and identify the topology IDs. Then, the 3D subgraph extraction is performed based on both 2D and 3D structure information. 

## Search for structure and substructure similarities in RAG-3D

RAG-3D's search engine for graph similarity considers two 2D graphs similar if they share the same pattern of connectivity among the vertices (same number of vertices and same $\lambda_2$). On the 3D level, the structural comparison is performed by a structural alignment of each possible 3D subgraph of a query structure to each of the 3D-graphs of the same topology ID. Structural differences between the optimally aligned graphs are measured by the Root Mean Square Deviation (RMSD) between the aligned vertices positions of the graphs. We rank the 3D-graphs based on their RMSDs to the query structure. The graph substructuring procedure is shown in ![Figure 3 in Zahran et al. NAR 2015](https://academic.oup.com/view-large/figure/81899335/gkv823fig3.jpeg).



## RAG-3D applications

Our 3D graphs, used here to search the RAG-3D database for molecules with the same topological characteristics, contain the tertiary information crucial for the search of functional similarity, and the search may be more efficient in graph space. Our graph comparison procedure based on 2D and 3D searches may thus aid the search for functional similarity. 

## Example
The python script **FindBestMatches.py** reads RNA secondary structure, 3D tree graph, loops and vertices information as input and obtains subgraphs and all-atom fragments that match with each subgraph. The input and output files of the 1DK1 example locate in ./Example/. The following command is used to run the script:

```
python FindBestMatches.py -query 1DK1 -input ./Examples/Input/ -output ./Examples/Output/ -database /path to RAG-3D-Database/ -graph graphID -matchVertexType True/False
```

Loops, Vertices, VertexTypes, Graph, BPSEQ files are needed in the input directory. Under the output directory 3 subfolders are created to store the results: Results, Subgraphs, and Junctions. The script also generates 2 subfolders: PML and MATCHES to store pymol script files and results of alighments, and the best matching results sorted with RMSD values. Best matches are searched against the RAG-3D database. The path to the database should be specified following **-database**. 

## References

Zahran, M., Sevim Bayrak C., Elmetwaly, S., Schlick, T. (2015) RAG-3D: a search tool for RNA 3D substructures. *Nucleic Acids Res.*, **19**, 9474-88.

Smit,S., Rother,K., Heringa,J. and Knight,R. (2008) From knotted to nested RNA structures: a variety of computational methods for pseudoknot removal. *RNA*, **14**, 410–416.

Laing,C., Jung,S., Kim,N., Elmetwaly,S., Zahran,M. and Schlick,T. (2013) Predicting helical topologies in RNA junctions as tree graphs. *PLoS One*, **8**, e71947.

Kim,N., Laing,C., Elmetwaly,S., Jung,S., Curuksu,J. and Schlick,T.
(2014) Graph-based sampling for approximating global helical topologies of RNA. *Proc. Natl. Acad. Sci. U.S.A.*, **111**, 4079–4084.

Kim,N., Zahran,M. and Schlick,T. (2015) Chapter five––computational prediction of riboswitch tertiary structures including pseudoknots by RAGTOP: a hierarchical graph sampling approach. *Method. Enzymol.*, **553**, 115–135.
