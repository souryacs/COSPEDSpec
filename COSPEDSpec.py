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
''' 

## Copyright 2014 Sourya Bhattacharyya and Jayanta Mukherjee.
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
  
  parser.add_option("-l", "--level", \
			  action="store_true", \
			  dest="level_info", \
			  default=False, \
			  help="if TRUE, number of levels of ancestor / descendent are considered during couplet relations. Default FALSE.")  
			  
  parser.add_option("-q", "--queues", \
			  type="int", \
			  action="store", \
			  dest="no_of_queues", \
			  default=1, \
			  help="1 - only a single max priority queue is used for storing the score metrics \
			  2 - two separate queues are used to store the conflicting and non conflicting taxa pairs and corresponding score metrics (default)")
  
  parser.add_option("-b", "--binary", \
			  type="int", \
			  action="store", \
			  dest="binary_suptr", \
			  default=2, \
			  help="0 - no binary refinement will be performed. \
			  1 - binary refinement on the basis of level based scoring and clustering (less runtime) \
			  2 - binary refinement on the basis of MRCA based extra lineage computation (more running time)")
  
  parser.add_option("-d", "--dfsref", \
			  action="store_true", \
			  dest="dfs_parent_refine", \
			  default=True, \
			  help="if TRUE, Multiple parent problem (C2) is tackled - before appying DFS for arbitrary parenting information, \
			it selects the most probable parent candidate.")  
    			    
  opts, args = parser.parse_args()
  return opts, args
  
  
##-----------------------------------------------------
''' main function '''
def main():  
  opts, args = parse_options()
  
  ROOTED_TREE = False 
  PRESERVE_UNDERSCORE = True 
  COST_UPDATE_LATEST = False 
  SINGLE_EDGE_TYPE_CONN_PRIORITY = True 
  INPUT_FILE_FORMAT = opts.inp_file_format
  INPUT_FILENAME = opts.INP_FILENAME
  NO_OF_QUEUES = opts.no_of_queues
  INCLUDE_ALL_RELN_FOR_CONFLICT_TAXA_PAIR = False 
  BINARY_SUPERTREE_OPTION = opts.binary_suptr
  OUTPUT_FILENAME = opts.OUT_FILENAME
  LEVEL_INFO_CONSIDER = opts.level_info
  DFS_PARENT_REFINE = opts.dfs_parent_refine
  
  global Output_Text_File
  
  if (INPUT_FILENAME == ""):
    print '******** THERE IS NO INPUT FILE SPECIFIED - RETURN **********'
    return
  else:
    print 'input filename: ', INPUT_FILENAME
    
  # according to the location of input filename
  # adjust the locations of the output files as well
  k = INPUT_FILENAME.rfind("/")
  if (k == -1):
    dir_of_inp_file = './'
  else:
    dir_of_inp_file = INPUT_FILENAME[:(k+1)]
  if (DEBUG_LEVEL > 1):
    print 'dir_of_inp_file: ', dir_of_inp_file  
      
  if (OUTPUT_FILENAME == ""):
    dir_of_curr_exec = dir_of_inp_file + 'COSPEDSpec'
    
    # comment - sourya
    ## update dir_of_curr_exec according to the settigs of input parameters
    #if (COST_UPDATE_LATEST == True):
      #dir_of_curr_exec = dir_of_curr_exec + '_C_1'
    #else:
      #dir_of_curr_exec = dir_of_curr_exec + '_C_0'
    #dir_of_curr_exec = dir_of_curr_exec + '_Q_' + str(NO_OF_QUEUES)
    #if (BINARY_SUPERTREE_OPTION == True):
      #dir_of_curr_exec = dir_of_curr_exec + '_B_1'
    #else:
      #dir_of_curr_exec = dir_of_curr_exec + '_B_0'
    #if (LEVEL_INFO_CONSIDER == True):
      #dir_of_curr_exec = dir_of_curr_exec + '_L_1'
    #else:
      #dir_of_curr_exec = dir_of_curr_exec + '_L_0'
    #if (DFS_PARENT_REFINE == True):
      #dir_of_curr_exec = dir_of_curr_exec + '_D_1'
    #else:
      #dir_of_curr_exec = dir_of_curr_exec + '_D_0'
    # end comment - sourya
      
    # append the current output directory in the text file
    Output_Text_File = dir_of_curr_exec + '/' + 'COSPEDSpec_Complete_Desription.txt'
    
    # create the directory
    if (os.path.isdir(dir_of_curr_exec) == False):
      mkdr_cmd = 'mkdir ' + dir_of_curr_exec
      os.system(mkdr_cmd)               
    
    if (DEBUG_LEVEL > 1):
      print 'Output_Text_File: ', Output_Text_File      
  else:
    # when we have specified the output file then corresponding directory becomes the dir_of_curr_exec
    k1 = OUTPUT_FILENAME.rfind("/")
    if (k1 == -1):
      dir_of_curr_exec = './'
    else:
      dir_of_curr_exec = OUTPUT_FILENAME[:(k1+1)]
    Output_Text_File = OUTPUT_FILENAME + '_complete_text_description'   
    if (DEBUG_LEVEL > 1):
      print 'Output_Text_File: ', Output_Text_File
    
  fp = open(Output_Text_File, 'w')    
  
  # comment - sourya
  #fp.write('\n ================ status of options ================= (1 means ON)')
  #fp.write('\n ROOTED_TREE: ' + str(ROOTED_TREE))
  #fp.write('\n PRESERVE_UNDERSCORE: ' + str(PRESERVE_UNDERSCORE))
  #fp.write('\n COST_UPDATE_LATEST: ' + str(COST_UPDATE_LATEST))
  #fp.write('\n SINGLE_EDGE_TYPE_CONN_PRIORITY: ' + str(SINGLE_EDGE_TYPE_CONN_PRIORITY))
  #fp.write('\n NO_OF_QUEUES: ' + str(NO_OF_QUEUES))
  #fp.write('\n INCLUDE_ALL_RELN: ' + str(INCLUDE_ALL_RELN_FOR_CONFLICT_TAXA_PAIR))
  #fp.write('\n LEVEL CONSIDERATION FOR COUPLET SCORING: ' + str(LEVEL_INFO_CONSIDER))
  #fp.write('\n MULTIPLE PARENT PROBLEM C2 solve: ' + str(DFS_PARENT_REFINE))
  #fp.write('\n STRICTLY BINARY SUPERTREE CREATE: ' + str(BINARY_SUPERTREE_OPTION))  
  #fp.write('\n ===>>>  processing the file now ======== ')
  # end comment - sourya
      
  # note the program beginning time 
  start_timestamp = time.time()
        
  #-------------------------------------  
  """ read the source trees collection and store it in a tree collection structure
  individual elements of this collection is thus a source tree """
  Gene_TreeList = Read_Gene_Data_Treelist(ROOTED_TREE, PRESERVE_UNDERSCORE, INPUT_FILE_FORMAT, INPUT_FILENAME)  
  
  #-------------------------------------    
  global number_of_taxa
  global tree_count
  
  ''' from the input source trees, define the pairwise taxa ordering and store it in convenient structure 
  also note the number of taxa (total) '''
  for k in Gene_TreeList:
    # "Curr_tree" contains the current source tree
    tree_count = tree_count + 1
    Curr_tree = Tree(k)    
    # update the input taxa list by generating leaves of the current tree 
    taxa_labels_curr_tree = Curr_tree.infer_taxa().labels()
    if (DEBUG_LEVEL > 0):
      fp.write('\n Tree no : ' + str(tree_count) + ' no of leaf nodes: ' + str(len(taxa_labels_curr_tree)))
    if (DEBUG_LEVEL > 2):
      fp.write('\n taxa set belonging to current tree: ' + str(taxa_labels_curr_tree))
  
    for i in range(len(taxa_labels_curr_tree)):
      ''' we also define one node structure "Taxa_Info_Dict" marked by a taxa
      current taxa is the key 
      for each taxa, one such structure is defined '''
      key1 = taxa_labels_curr_tree[i]
      if key1 not in Taxa_Info_Dict:
	if (DEBUG_LEVEL > 2):
	  fp.write('\n new node in Taxa_Info_Dict : ' + str(key1))
	Taxa_Info_Dict.setdefault(key1, Single_Taxa())
	COMPLETE_INPUT_TAXA_LIST.append(taxa_labels_curr_tree[i])
	number_of_taxa = number_of_taxa + 1	# increment the counter
  
    ''' now for each pair of taxa, we define another structure containing their relation
    that is, they are ancestor / descendent, or no relation at all
    the structure is formed according to the input source trees
    so one pair of taxa may have multiple types of relation (according to different source trees)
    also any relation may have multiple support (no of source trees confirming the relation) '''
    for i in range(len(taxa_labels_curr_tree) - 1):
      for j in range(i+1, len(taxa_labels_curr_tree)):
	# keys are formed by leaf (taxa) pairs
	key1 = (taxa_labels_curr_tree[i], taxa_labels_curr_tree[j])
	key2 = (taxa_labels_curr_tree[j], taxa_labels_curr_tree[i])
	''' we define the relationship between these two taxa 
	basically we mark particular relation type '''  	
	# we check whether the pair is already existing in the list of taxa pairs
	if (key1 not in TaxaPair_Reln_Dict) and (key2 not in TaxaPair_Reln_Dict):
	  if (DEBUG_LEVEL > 2):
	    fp.write('\n TaxaPair_Reln_Dict has new element' + str(key1))
	  TaxaPair_Reln_Dict.setdefault(key1, Reln_TaxaPair(UNDEFINED_EDGE))
	  # this function defines the taxa pair relationship
	  # with respect to particular source tree
	  DefineLeafPairReln(taxa_labels_curr_tree[i], taxa_labels_curr_tree[j], Curr_tree, \
	    LEVEL_INFO_CONSIDER, COST_UPDATE_LATEST)
	elif (key1 in TaxaPair_Reln_Dict):
	  DefineLeafPairReln(taxa_labels_curr_tree[i], taxa_labels_curr_tree[j], Curr_tree, \
	    LEVEL_INFO_CONSIDER, COST_UPDATE_LATEST)
	else:
	  DefineLeafPairReln(taxa_labels_curr_tree[j], taxa_labels_curr_tree[i], Curr_tree, \
	    LEVEL_INFO_CONSIDER, COST_UPDATE_LATEST)
  
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
  
  fp.close()
  
  # we print the original connection status for all the tree nodes
  if 0:	#(DEBUG_LEVEL > 2):
    for l in Taxa_Info_Dict:
      #print 'printing information for the Taxa_Info_Dict key ', l
      Taxa_Info_Dict[l]._PrintOriginalTaxaInfo(l, Output_Text_File)      
        
  data_read_timestamp = time.time()	# note the timestamp
    
  #------------------------------------------------------------
  """ this section contains the code to estimate the supertree 
  from the input source tree list information """
  
  ''' we also calculate the connection value between each pair of nodes in the output tree
  the value defines the majority of the edge type that is between those 2 nodes '''
  for l in TaxaPair_Reln_Dict:
    ''' calculate the cost associated with each node pair connection for different edge types 
    single_edge_exist means that only one type of edge connection is set between these two nodes '''
    if (SINGLE_EDGE_TYPE_CONN_PRIORITY):
      # detect if only one type of connection exists, during setting the priority values of different edge connections
      single_edge_exist_list = TaxaPair_Reln_Dict[l]._SetConnPrVal(True)
      single_edge_exist = single_edge_exist_list[0]
      edge_type = single_edge_exist_list[1]
    else:
      # detect if only one type of connection exists, during setting the priority values of different edge connections
      single_edge_exist = TaxaPair_Reln_Dict[l]._SetConnPrVal(False)
      edge_type = NO_EDGE	# only used when single_edge_exist is 1

    """ we calculate the cost metric value between individual pairs of taxa
    previously the cost metric was equal to the priority metric
    now we change it to make it a product of the frequency and the priority metrics """
    TaxaPair_Reln_Dict[l]._SetCostMetric()
  
    #------------------------------------------------------------
    """ also update the cost value for individual elements in the 
    list "Cost_List_Taxa_Pair_Multi_Reln or Cost_List_Taxa_Pair_Single_Reln"
    each list element contains the following values:
    1) taxa1 and taxa2    
    3) edge type (one at a time - so there will be 4 entries for each node pair)
    4) edge cost (the cost associated with one particular edge type """
    
    #-------------------
    # for taxa pairs supported by only one source tree, if -a option is false
    # then only one relation is included
    # otherwise all the relations are included
    #-------------------
    
    if (single_edge_exist == 0):
      for edge_type in range(4):
	if (INCLUDE_ALL_RELN_FOR_CONFLICT_TAXA_PAIR == True):
	  # for the conflicting taxa pair (multiple relations between the taxa pair are supported by the source trees)
	  # there all the relations were included in the score metric queue
	  sublist = [l[0], l[1], edge_type, TaxaPair_Reln_Dict[l]._GetEdgeCost_ConnReln(edge_type)]
	  Cost_List_Taxa_Pair_Multi_Reln.append(sublist)
	else:
	  # now we add only those relations between the taxa pair which are supported by at least one source tree
	  if (TaxaPair_Reln_Dict[l]._GetEdgeWeight(edge_type) > 0):
	    sublist = [l[0], l[1], edge_type, TaxaPair_Reln_Dict[l]._GetEdgeCost_ConnReln(edge_type)]
	    Cost_List_Taxa_Pair_Multi_Reln.append(sublist)
    else:
      # for single edge type, assign that particular connection
      sublist = [l[0], l[1], edge_type, TaxaPair_Reln_Dict[l]._GetEdgeCost_ConnReln(edge_type)]
      """ this connection is only possible between the current examined taxa pairs
      if SINGLE_EDGE_TYPE_CONN_PRIORITY is TRUE then the information is placed in Cost_List_Taxa_Pair_Single_Reln
      otherwise it is placed in Cost_List_Taxa_Pair_Multi_Reln """
      if (SINGLE_EDGE_TYPE_CONN_PRIORITY) and (NO_OF_QUEUES == 2):
	Cost_List_Taxa_Pair_Single_Reln.append(sublist)
      else:
	Cost_List_Taxa_Pair_Multi_Reln.append(sublist)
    
  #------------------------------------------------------------  
  if (DEBUG_LEVEL >= 2):
    for l in TaxaPair_Reln_Dict:
      #print 'printing info for the TaxaPair_Reln_Dict key: ', l
      TaxaPair_Reln_Dict[l]._PrintRelnInfo(l, Output_Text_File)
        
  """ here we allocate the list of clusters
  initially all the clusters contain one taxa
  each of the cluster has the index of the corresponding taxa in the COMPLETE_INPUT_TAXA_LIST """
  for i in range(len(COMPLETE_INPUT_TAXA_LIST)):
    Create_Cluster_Taxa_Label(i, COMPLETE_INPUT_TAXA_LIST[i])
  
  #------------------------------------------------------------
  """ we initialize the Reachability_Graph_Mat
  for all the clusters of nodes possible, this indicate the edges between a cluster pair
  initially the number of clusters are thought of equal as the number of taxa
  gradually, as the taxa subsets are merged in a cluster, corresponding count is decreased 
  convention -  out edge from the first node to the second node 
  this is a numpy 2D array 
  values Mat[x][y] = 1 means x->y
  Mat[x][y] = Mat[y][x] = 2 means x and y are connected via NO EDGE """
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
  """ now we have to sort the Cost_List_Node_Pair according to the edge cost value 
  that is the 4th field of the sublist 
  we use custom sorting operation """
  Sort_Cost_List_Initial(Cost_List_Taxa_Pair_Multi_Reln)
  
  # if there is provision to include single connectivity edges then we sort that list also
  if (SINGLE_EDGE_TYPE_CONN_PRIORITY) and (NO_OF_QUEUES == 2):
    Sort_Cost_List_Initial(Cost_List_Taxa_Pair_Single_Reln)
  
  data_initialize_timestamp = time.time()	# note the timestamp
  
  #------------------------------------------------------------
  """ if we have stored taxa pairs depicting single relation instance 
  then we first process that corresponding queue """
  if (SINGLE_EDGE_TYPE_CONN_PRIORITY) and (NO_OF_QUEUES == 2):
    Reachability_Graph_Mat = Proc_Queue(Reachability_Graph_Mat, COST_UPDATE_LATEST, 1, Output_Text_File, BINARY_SUPERTREE_OPTION, DFS_PARENT_REFINE)
  
  # then we process the queue containing the taxa pairs depicting multiple relation instance 
  Reachability_Graph_Mat = Proc_Queue(Reachability_Graph_Mat, COST_UPDATE_LATEST, 0, Output_Text_File, BINARY_SUPERTREE_OPTION, DFS_PARENT_REFINE)
  
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
  """ now perform the transitive reduction of the closure formed by connection of the cluster of nodes in the above operation
  this is required to handle the following scenario:
  suppose, there exists a case such that A->C, B->C and A->B
  then in the final graph, only A->B and B->C information needs to be preserved
  in order to form the DAG 
  
  this is the solution of the problem C1 mentioned in the supertree manuscript"""
  CompressDirectedGraph(Reachability_Graph_Mat)
  
  #------------------------------------------------------------
  #""" now instead of arbitrary assignment of the parent node for individual clusters 
  #we assign parent node according to the source tree relationships
  #this will solve the multiple parent problem C2 as discussed in the manuscript 
  
  #this is a new addition and marked under the DFS based parent refinement option """
  if (DFS_PARENT_REFINE == True):
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
  ''' now this section constructs the supertree from the generated DAG 
  this is performed by repeatedly extracting the nodes with minimum indegree
  basically we first form a string which represents the supertree '''
  no_of_components = 0	# for forest
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
  
  # with the final tree string, finally generate the tree result 
  Final_Supertree_Str = '(' + Final_Supertree_Str + ')'
  
  #--------------------------------------------------------------
  # add - sourya
  fp = open(Output_Text_File, 'a')
  fp.write('\n --- original supertree as newick string --- ' + Final_Supertree_Str) 
  Final_Supertree_Str = Remove_Extra_Paranthesis(Final_Supertree_Str)  
  fp.write('\n --- after removing extra paranthesis -- supertree as newick string --- ' + Final_Supertree_Str) 
  # end add - sourya
  #--------------------------------------------------------------
  # final timestamp
  data_process_timestamp = time.time()      
  
  # read the supertree (without branch length information)
  Supertree_without_branch_len = dendropy.Tree.get_from_string(Final_Supertree_Str, schema="newick", \
								preserve_underscores=PRESERVE_UNDERSCORE, \
								default_as_rooted=True)          
  
  # add - sourya
  Supertree_without_branch_len.update_splits(delete_outdegree_one=True)
  fp.write('\n --- after update splits -- supertree : ' + Supertree_without_branch_len.as_newick_string())
  fp.close()
  # end add - sourya
  
  #--------------------------------------------------------------
  # add - sourya  
  if (BINARY_SUPERTREE_OPTION > 0):
    fp = open(Output_Text_File, 'a')
    if (BINARY_SUPERTREE_OPTION == 1):
      fp.write('\n --- user provided option for producing strict binary supertree - refinement with level based scoring computation')
    else:
      fp.write('\n --- user provided option for producing strict binary supertree - refinement with MRCA based extra lineage computation')
    fp.close()
    # this function removes all multifurcating clusters and produces binary tree (except problem C3)
    Refine_Supertree_Binary_Form(Gene_TreeList, Supertree_without_branch_len, Output_Text_File, BINARY_SUPERTREE_OPTION)
    fp = open(Output_Text_File, 'a')
    fp.write('\n --- after binary refinement --- output tree without branch length (in newick format): ' + Supertree_without_branch_len.as_newick_string())    
    fp.close()
  else:
    fp = open(Output_Text_File, 'a')
    fp.write('\n --- user did not provide option for producing strict binary supertree - so output tree can be non-binary')
    fp.close()
  
  # final timestamp
  binary_refinement_timestamp = time.time()      
  # end add - sourya  
  
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
  if (len(EDGE_PROCESSED_LIST) > 0):
    EDGE_PROCESSED_LIST[:] = []
  if (len(COMPLETE_INPUT_TAXA_LIST) > 0):
    COMPLETE_INPUT_TAXA_LIST[:] = []
  if (len(CURRENT_CLUST_IDX_LIST) > 0):
    CURRENT_CLUST_IDX_LIST[:] = []
  
  # free the reachability graph (numpy array)
  del Reachability_Graph_Mat
    
#-----------------------------------------------------

if __name__ == "__main__":
    main() 
  
