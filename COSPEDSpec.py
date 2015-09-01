#!/usr/bin/env python

##---------------------------------------------
''' 
this program is used to generate a spcies tree from a set of gene trees
the gene trees generally associate conflicts among the constituent genes (representing individual taxa)
in terms of the topological relationships
species tree follows supertree construction approach to produce species tree with high performance index

Author: Sourya Bhattacharyya
Dept of CSE, IIT Kharagpur
V1.0 - 07.04.2014 - basic code
V2.0 - 01.05.2014 - modified the clustering and conflict detection routine
V3.0 - 05.09.2014 - incorporated binary supertree (fully resolved)
V4.0 - 07.11.2014 - storage space clean, cost update stop, code clean
V5.0 - 14.04.2015 - 1) phylogenetic header library use for coding, 2) optimization in tree reading procedures
V6.0 - 01.09.2015 - code clean and comments
''' 

## Copyright 2014, 2015 Sourya Bhattacharyya and Jayanta Mukherjee.
## All rights reserved.
##
## See "LICENSE.txt" for terms and conditions of usage.
##
##---------------------------------------------

import Header
from Header import *
import Cost_Update
from Cost_Update import *
import Process_Queues
from Process_Queues import *
import ReachGraph_Update
from ReachGraph_Update import *
import UtilFunc
from UtilFunc import *
import RefineSpeciesTree
from RefineSpeciesTree import *

##-----------------------------------------------------
# this function is useful to parse various options for input data processing
def parse_options():  
  parser = OptionParser()
    
  parser.add_option("-I", "--INPFILE", \
			  type="string", \
			  action="store", \
			  dest="INP_FILENAME", \
			  default="", \
			  help="name of the input file containing gene trees")
			  
  parser.add_option("-O", "--OUTFILE", \
			  type="string", \
			  action="store", \
			  dest="OUT_FILENAME", \
			  default="", \
			  help="name of the output file to contain target species tree")  
  
  parser.add_option("-p", "--inpform", \
			  type="int", \
			  action="store", \
			  dest="inp_file_format", \
			  default=1, \
			  help="1 - input file format is NEWICK (default) \
			  2 - input file format is NEXUS")
  			  
  #parser.add_option("-q", "--queues", \
			  #type="int", \
			  #action="store", \
			  #dest="no_of_queues", \
			  #default=2, \
			  #help="1 - only a single max priority queue is used for storing the score metrics \
			  #2 - two separate queues are used to store the conflicting and non conflicting taxa pairs and corresponding score metrics (default)")
   
  opts, args = parser.parse_args()
  return opts, args
  
  
##-----------------------------------------------------
# main function
def main():  
  opts, args = parse_options()
  
  ROOTED_TREE = False 
  PRESERVE_UNDERSCORE = True  
  if (opts.inp_file_format == 1):
    INPUT_FILE_FORMAT = 'newick'
  else:
    INPUT_FILE_FORMAT = 'nexus'
  INPUT_FILENAME = opts.INP_FILENAME
  NO_OF_QUEUES = 2	#opts.no_of_queues  
  OUTPUT_FILENAME = opts.OUT_FILENAME
  
  global Output_Text_File
  
  if (INPUT_FILENAME == ""):
    print '******** THERE IS NO INPUT FILE SPECIFIED - RETURN **********'
    return
  else:
    print 'input filename: ', INPUT_FILENAME
    
  """ 
  according to the location of input filename
  adjust the locations of the output files as well
  """
  k = INPUT_FILENAME.rfind("/")
  if (k == -1):
    dir_of_inp_file = './'
  else:
    dir_of_inp_file = INPUT_FILENAME[:(k+1)]
  if (DEBUG_LEVEL > 1):
    print 'dir_of_inp_file: ', dir_of_inp_file  
      
  if (OUTPUT_FILENAME == ""):
    dir_of_curr_exec = dir_of_inp_file + 'COSPEDSpec'
    """ 
    create the directory
    """
    if (os.path.isdir(dir_of_curr_exec) == False):
      mkdr_cmd = 'mkdir ' + dir_of_curr_exec
      os.system(mkdr_cmd)                     
    """ 
    location of the output text file (where results would be written)
    """
    Output_Text_File = dir_of_curr_exec + '/' + 'COSPEDSpec_Complete_Desription.txt'
    print 'Output_Text_File: ', Output_Text_File      
  else:
    """ 
    when we have specified the output file 
    then corresponding directory becomes the dir_of_curr_exec
    """
    k1 = OUTPUT_FILENAME.rfind("/")
    if (k1 == -1):
      dir_of_curr_exec = './'
    else:
      dir_of_curr_exec = OUTPUT_FILENAME[:(k1+1)]
    Output_Text_File = OUTPUT_FILENAME + '_complete_text_description'   
    print 'Output_Text_File: ', Output_Text_File
    
  # open the output text file
  fp = open(Output_Text_File, 'w')    
  
  #fp.write('\n ================ status of options ================= (1 means ON)')
  #fp.write('\n ROOTED_TREE: ' + str(ROOTED_TREE))
  #fp.write('\n PRESERVE_UNDERSCORE: ' + str(PRESERVE_UNDERSCORE))
  fp.write('\n NO_OF_QUEUES: ' + str(NO_OF_QUEUES))
  fp.write('\n ===>>>  processing the file now ======== ')
      
  # note the program beginning time 
  start_timestamp = time.time()
        
  #-------------------------------------  
  """ 
  read the source trees collection and store it in a treelist
  individual elements of this collection is a source tree 
  """
  Gene_TreeList = Read_Input_Treelist(ROOTED_TREE, PRESERVE_UNDERSCORE, INPUT_FILE_FORMAT, INPUT_FILENAME)  
  
  #-------------------------------------
  tree_count = 0
  
  """ 
  from the input source trees, note the number of taxa (total)
  and also define the class instances corresponding to single taxa
  """
  for tr_idx in range(len(Gene_TreeList)):
    tree_count = tree_count + 1
    taxa_labels_curr_tree = Gene_TreeList[tr_idx].infer_taxa().labels()
    if (DEBUG_LEVEL > 1):
      fp.write('\n Tree no : ' + str(tree_count) +  'no of leaf nodes: ' + str(len(taxa_labels_curr_tree)))
    if (DEBUG_LEVEL > 2):
      fp.write('\n taxa set belonging to current tree: ' + str(taxa_labels_curr_tree))
    for i in range(len(taxa_labels_curr_tree)):
      if taxa_labels_curr_tree[i] not in COMPLETE_INPUT_TAXA_LIST:
	COMPLETE_INPUT_TAXA_LIST.append(taxa_labels_curr_tree[i])
  
  # we also define one structure "Taxa_Info_Dict" marked by a taxa
  for label in COMPLETE_INPUT_TAXA_LIST:
    Taxa_Info_Dict.setdefault(label, Single_Taxa())
  
  # now process individual trees to find the couplet relations of those trees
  for tr_idx in range(len(Gene_TreeList)):
    DeriveCoupletRelations(Gene_TreeList[tr_idx])
  
  number_of_taxa = len(COMPLETE_INPUT_TAXA_LIST)
  
  if (DEBUG_LEVEL >= 0):
    fp.write('\n  total no of taxa: ' + str(number_of_taxa))
  if (DEBUG_LEVEL > 2):
    fp.write('\n len Taxa_Info_Dict: ' + str(len(Taxa_Info_Dict)))
    fp.write('\n len COMPLETE_INPUT_TAXA_LIST: ' + str(len(COMPLETE_INPUT_TAXA_LIST))) 
    fp.write('\n len TaxaPair_Reln_Dict : ' + str(len(TaxaPair_Reln_Dict)))
    for i in range(len(COMPLETE_INPUT_TAXA_LIST)):
      fp.write('\n COMPLETE_INPUT_TAXA_LIST idx : ' + str(i) + ' element: ' + str(COMPLETE_INPUT_TAXA_LIST[i]))
  
  if (DEBUG_LEVEL > 2):
    for l in Taxa_Info_Dict:
      fp.write('\n taxa index in COMPLETE_INPUT_TAXA_LIST: ' + str(COMPLETE_INPUT_TAXA_LIST.index(l)))
  
  # close the output text file
  fp.close()
  
  # we print the original connection status for all the tree nodes
  if 0:	#(DEBUG_LEVEL > 2):
    for l in Taxa_Info_Dict:
      #print 'printing information for the Taxa_Info_Dict key ', l
      Taxa_Info_Dict[l]._PrintOriginalTaxaInfo(l, Output_Text_File)      
        
  data_read_timestamp = time.time()	# note the timestamp
    
  #------------------------------------------------------------
  """ 
  this section contains the code to estimate the supertree 
  from the input treelist 
  """  
  for l in TaxaPair_Reln_Dict:
    """ 
    calculate the support score and priority measures for individual couplets
    and for individual relations
    """
    """
    single_edge_exist: if TRUE, means that only one type of relation is supported (with respect to input trees) between this couplet
    detection of it during setting the priority values of different relations
    """
    single_edge_exist_list = TaxaPair_Reln_Dict[l]._SetConnPrVal(True)
    single_edge_exist = single_edge_exist_list[0]
    edge_type = single_edge_exist_list[1]

    """ 
    we calculate the support score between individual couplet
    it is the product of frequency and the priority measures 
    """
    TaxaPair_Reln_Dict[l]._SetCostMetric()
  
    #------------------------------------------------------------
    """ 
    support score and corresponding couplets will be placed in one of the following priority queues
    1) Cost_List_Taxa_Pair_Multi_Reln: for conflicting couplets (more than one relations are supported)
    2) Cost_List_Taxa_Pair_Single_Reln: for non-conflicting couplets + when we use two different queues 
    each element of the queue contains the following fields:
    1) taxa1 and taxa2    
    2) relation type (r1 to r4)
    3) support score for that relation 
    """    
    if (single_edge_exist == 0):
      for edge_type in range(4):
	""" 
	note: we add only those relations (between a taxa pair) in a queue 
	which are supported by at least one source tree
	that is, their frequency is non zero
	"""
	if (TaxaPair_Reln_Dict[l]._GetEdgeWeight(edge_type) > 0):
	  sublist = [l[0], l[1], edge_type, TaxaPair_Reln_Dict[l]._GetEdgeCost_ConnReln(edge_type)]
	  Cost_List_Taxa_Pair_Multi_Reln.append(sublist)
    else:
      """ 
      if the current relation is the strict consensus and only supported relation for this couplet, 
      copy the relation information in appropriate queue
      """
      sublist = [l[0], l[1], edge_type, TaxaPair_Reln_Dict[l]._GetEdgeCost_ConnReln(edge_type)]
      """ 
      if we have provision for using two separate queues, then we place this relation in Cost_List_Taxa_Pair_Single_Reln
      else we place this relation in Cost_List_Taxa_Pair_Multi_Reln
      """
      if (NO_OF_QUEUES == 2):
	Cost_List_Taxa_Pair_Single_Reln.append(sublist)
      else:
	Cost_List_Taxa_Pair_Multi_Reln.append(sublist)
    
  #------------------------------------------------------------  
  if (DEBUG_LEVEL >= 2):
    for l in TaxaPair_Reln_Dict:
      #print 'printing info for the TaxaPair_Reln_Dict key: ', l
      TaxaPair_Reln_Dict[l]._PrintRelnInfo(l, Output_Text_File)
        
  """ 
  Here we allocate the list / storage space of taxa clusters
  each taxa cluster is supposed to store the taxa subsets related via relation r3 (simultaneous speciation)
  initially all the clusters contain one taxa
  each of the cluster has the index of the corresponding taxa in the COMPLETE_INPUT_TAXA_LIST 
  """
  for i in range(len(COMPLETE_INPUT_TAXA_LIST)):
    Create_Cluster_Taxa_Label(i, COMPLETE_INPUT_TAXA_LIST[i])
  
  #------------------------------------------------------------
  """ 
  we initialize the Reachability_Graph_Mat
  dimension: N X N where N = number of taxa clusters = number of taxa (initially)
  this is a numpy 2D array 
  values Mat[x][y] = 1 means x->y
  Mat[x][y] = Mat[y][x] = 2 means x and y are connected via relation r4
  as taxa subsets are clustered, no of distinct taxa clusters reduces
  so, this matrix dimension is also decreased
  """
  fp = open(Output_Text_File, 'a')
  
  Reachability_Graph_Mat = numpy.zeros((number_of_taxa, number_of_taxa), dtype=numpy.int)
  if (DEBUG_LEVEL > 0):
    fp.write('\n shape of Reachability_Graph_Mat: ' + str(numpy.shape(Reachability_Graph_Mat)))
  
  # this information is printed to know the maximum possible iterations that the while loops will undergo
  if (DEBUG_LEVEL > 1):
    fp.write('\n =========== max connection pair ============= : ' + str((len(Cost_List_Taxa_Pair_Single_Reln) + len(Cost_List_Taxa_Pair_Multi_Reln))))      
    fp.write('\n len Cost_List_Taxa_Pair_Single_Reln: ' + str(len(Cost_List_Taxa_Pair_Single_Reln)))
    fp.write('\n len Cost_List_Taxa_Pair_Multi_Reln : ' + str(len(Cost_List_Taxa_Pair_Multi_Reln)))
  
  fp.close()
  #------------------------------------------------------------
  """ 
  sort the priority queues according to the support score value of individual relations
  4th field of individual sublist denotes the support score
  we use max priority queue implementation for storing and sorting these queues
  """
  # sorting the conflicting queue
  Sort_Cost_List_Initial(Cost_List_Taxa_Pair_Multi_Reln)
  # if we are allowed to use two distinct queues, then sort the non conflicting queue also
  if (NO_OF_QUEUES == 2):
    Sort_Cost_List_Initial(Cost_List_Taxa_Pair_Single_Reln)
  
  data_initialize_timestamp = time.time()	# note the timestamp
  
  #------------------------------------------------------------
  """
  process the non-conflicting queue first, provided we use this queue
  """
  if (NO_OF_QUEUES == 2):
    Reachability_Graph_Mat = Proc_Queue(Reachability_Graph_Mat, 1, Output_Text_File)
  
  """ 
  then we process the conflicting queue (couplets having multiple relations supported)
  """
  Reachability_Graph_Mat = Proc_Queue(Reachability_Graph_Mat, 0, Output_Text_File)
  
  # we print the final connection status for all the tree nodes
  if (DEBUG_LEVEL > 2):
    for l in Taxa_Info_Dict:
      #print 'printing information for the Taxa ', l
      Taxa_Info_Dict[l]._PrintFinalTaxaInfo(l, Output_Text_File) 
  
  # print the cluster information 
  if (DEBUG_LEVEL > 2):
    fp = open(Output_Text_File, 'a')
    fp.write('\n **** total number of clusters: ' + str(len(CURRENT_CLUST_IDX_LIST)))
    fp.write('\n CURRENT_CLUST_IDX_LIST contents: ')
    fp.write(str(CURRENT_CLUST_IDX_LIST))
    fp.write('\n ========== cluster information after reachability graph generation =============')
    fp.close()
    for i in Cluster_Info_Dict:
      #print 'printing the information for cluster node: ', i
      Cluster_Info_Dict[i]._PrintClusterInfo(i, Output_Text_File)
  
  # note the timestamp
  reachability_graph_form_timestamp = time.time()  
  #------------------------------------------------------------
  """ 
  after processing above queue, taxa subsets (clusters) are formed (r3 relation)
  and their mutual connectivity information (relations r1, r2 and r4) are established
  now perform the transitive reduction of taxa clusters
  to handle the following scenario:
  suppose, there exists a case such that A->C, B->C and A->B
  then in the final graph, only A->B and B->C information needs to be preserved
  in order to form the DAG 
  *** Note: this is the solution of the problem C1 (mentioned in the COSPEDTree manuscript)
  """
  CompressDirectedGraph(Reachability_Graph_Mat)
  
  #------------------------------------------------------------
  """ 
  *** addition: new to COSPEDSpec ****
  instead of arbitrary assignment of the parent node for individual clusters 
  we assign parent node according to the source tree relationships
  this will solve the multiple parent problem C2 (as discussed in the COSPEDSpec manuscript) 
  """
  SolveMultipleParentC2Problem(Output_Text_File)
  
  #------------------------------------------------------------
  # print the cluster information 
  if (DEBUG_LEVEL > 2):
    fp = open(Output_Text_File, 'a')
    fp.write('\n **** total number of clusters: ' + str(len(CURRENT_CLUST_IDX_LIST)))
    fp.write('\n CURRENT_CLUST_IDX_LIST contents: ')
    fp.write(str(CURRENT_CLUST_IDX_LIST))    
    fp.write('\n ========== cluster information after transitive reduction =============')
    fp.close()
    for i in Cluster_Info_Dict:
      #print 'printing the information for cluster node: ', i
      Cluster_Info_Dict[i]._PrintClusterInfo(i, Output_Text_File)
      
  # note the timestamp
  cluster_of_node_refine_species_timestamp1 = time.time()  
  #----------------------------------------------  
  """ 
  construct supertree from the generated DAG 
  scheme: repeatedly extract the nodes (taxa clusters) with minimum indegree
  print taxa information constituent within this cluster
  maintain appropriate paranthesis according to the newick format
  """
  """ 
  the variable no_of_components is used in association with the used depth first technique 
  if no_of_components > 1, it signifies that a forest has been created
  in other words, the problem C3 (no parent problem exists)
  """
  no_of_components = 0	
  while (1):
    root_clust_node_idx = Extract_Node_Min_Indeg(len(CURRENT_CLUST_IDX_LIST))
    if (root_clust_node_idx == -1):
      break
    Tree_Str = PrintNewick(root_clust_node_idx)
    no_of_components = no_of_components + 1
    if (no_of_components == 1):	# first component
      Final_Supertree_Str = Tree_Str
    else:
      Final_Supertree_Str = Final_Supertree_Str + ',' + Tree_Str
  
  """ 
  with the final tree string, finally generate the tree result 
  this is also required to tackle the creation of a forest during the DFS (no parent problem)
  this procedure adheres to the basic COSPEDTree mechanism
  """
  Final_Supertree_Str = '(' + Final_Supertree_Str + ')'
  
  #--------------------------------------------------------------
  """
  Here, we first modify the generated newick string (representing the supertree)
  to delete unnecessary internal edges, etc.
  *** note: we did not use dendropy function since we encountered some errors
  """
  fp = open(Output_Text_File, 'a')
  fp.write('\n --- original supertree as newick string --- ' + Final_Supertree_Str) 
  Final_Supertree_Str = Remove_Extra_Paranthesis(Final_Supertree_Str)  
  fp.write('\n --- after removing extra paranthesis -- supertree as newick string --- ' + Final_Supertree_Str) 
  
  # final timestamp
  data_process_timestamp = time.time()      
  #--------------------------------------------------------------
  
  # read the supertree (without branch length information)
  Supertree_without_branch_len = dendropy.Tree.get_from_string(Final_Supertree_Str, schema="newick", \
								preserve_underscores=PRESERVE_UNDERSCORE, \
								default_as_rooted=True)          
  
  ## add - sourya
  #Supertree_without_branch_len.update_splits(delete_outdegree_one=True)
  #fp.write('\n --- after update splits -- supertree : ' + Supertree_without_branch_len.as_newick_string())
  #fp.close()
  ## end add - sourya
  
  #--------------------------------------------------------------
  fp = open(Output_Text_File, 'a')
  fp.write('\n --- user provided option for producing strict binary supertree')
  fp.close()
  
  """ 
  this function removes all multifurcating clusters and produces binary tree 
  it also solves the problem C3, as mentioned in the manuscript
  """
  Refine_Supertree_Binary_Form(Gene_TreeList, Supertree_without_branch_len, Output_Text_File)
  
  fp = open(Output_Text_File, 'a')
  fp.write('\n --- after binary refinement --- output tree without branch length (in newick format): ' + Supertree_without_branch_len.as_newick_string())    
  fp.close()
  
  # final timestamp
  binary_refinement_timestamp = time.time()      
  
  #--------------------------------------------------------------  
  # write this tree on a separate text file
  if (OUTPUT_FILENAME == ""):
    out_treefilename = dir_of_curr_exec + '/' + 'outtree_Newick.tre'
  else:
    out_treefilename = OUTPUT_FILENAME
  
  # we write the unweighted supertree
  outfile = open(out_treefilename, 'w')
  outfile.write(Supertree_without_branch_len.as_newick_string())
  outfile.close()
  
  # we write the time associated with the execution of this method
  fp = open(Output_Text_File, 'a')
  fp.write('\n \n\n ===============>>>>>>>>>>>>>>> TIME COMPLEXITY OF THE METHOD (in seconds) ')
  fp.write('\n \n reading the data: ' + str(data_read_timestamp - start_timestamp) + \
	'\n initialization of the structure: ' + str(data_initialize_timestamp - data_read_timestamp) + \
	'\n formation of the reachability graph (cluster) (after loop): ' + \
	      str(reachability_graph_form_timestamp - data_initialize_timestamp) + \
	'\n multiple parent (related) problem: ' + \
	      str(cluster_of_node_refine_species_timestamp1 - reachability_graph_form_timestamp) + \
	'\n newick string formation for non branch length trees: ' + \
	      str(data_process_timestamp - cluster_of_node_refine_species_timestamp1) + \
	'\n refining tree for producing strict binary supertree: ' + \
	      str(binary_refinement_timestamp - data_process_timestamp))
  fp.write('\n \n Total time taken (in seconds) : ' + str(binary_refinement_timestamp - start_timestamp))
  fp.close()
  
  #--------------------------------------------------------------  
  # delete the storage variables associated with the current execution 
  
  # clear the dictionaries
  Cluster_Info_Dict.clear()
  Taxa_Info_Dict.clear()
  TaxaPair_Reln_Dict.clear()
  
  # clear the lists associated
  if (len(Cost_List_Taxa_Pair_Multi_Reln) > 0):
    Cost_List_Taxa_Pair_Multi_Reln[:] = []
  if (len(Cost_List_Taxa_Pair_Single_Reln) > 0):
    Cost_List_Taxa_Pair_Single_Reln[:] = []
  if (len(COMPLETE_INPUT_TAXA_LIST) > 0):
    COMPLETE_INPUT_TAXA_LIST[:] = []
  if (len(CURRENT_CLUST_IDX_LIST) > 0):
    CURRENT_CLUST_IDX_LIST[:] = []
  
  # free the reachability graph (numpy array)
  del Reachability_Graph_Mat
    
#-----------------------------------------------------
if __name__ == "__main__":
    main() 
  
