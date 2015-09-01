#!/usr/bin/env python

import Header
from Header import *

##-----------------------------------------------------
# this function computes the score (ancestor relation) from clust1 to clust2
def ComputeScore(clust1, clust2, Output_Text_File):
  
  if (DEBUG_LEVEL > 2):
    fp = open(Output_Text_File, 'a')
    fp.write('\n ===>>> Score --- cluster 1 - taxa set: ' + str(Cluster_Info_Dict[clust1]._GetSpeciesList())\
      + ' cluster 2 taxa set ' + str(Cluster_Info_Dict[clust2]._GetSpeciesList()))
    
  score_val = 0
  for t1 in Cluster_Info_Dict[clust1]._GetSpeciesList():
    for t2 in Cluster_Info_Dict[clust2]._GetSpeciesList():
      key1 = (t1, t2)
      key2 = (t2, t1)
      if key1 in TaxaPair_Reln_Dict:
	# comment - sourya
	# previously the priority metric based scoring was employed to determine the scoring 
	# among two taxa clusters
	score_val = score_val + TaxaPair_Reln_Dict[key1]._GetConnPrVal(DIRECTED_OUT_EDGE)
	# add - sourya
	# now we employ frequency based scoring mechanism, to satisfy maximum agreement property
	#score_val = score_val + TaxaPair_Reln_Dict[key1]._GetEdgeWeight(DIRECTED_OUT_EDGE)
	# end add - sourya
      elif key2 in TaxaPair_Reln_Dict:
	# comment - sourya
	# previously the priority metric based scoring was employed to determine the scoring 
	# among two taxa clusters	
	score_val = score_val + TaxaPair_Reln_Dict[key2]._GetConnPrVal(DIRECTED_IN_EDGE)
	# add - sourya
	# now we employ frequency based scoring mechanism, to satisfy maximum agreement property
	#score_val = score_val + TaxaPair_Reln_Dict[key2]._GetEdgeWeight(DIRECTED_IN_EDGE)
	# end add - sourya
      else:
	if (DEBUG_LEVEL > 2):
	  fp.write('\n score compute -- key pair ' + str(t1) + ',' + str(t2) + ' does not exist ')

  if (DEBUG_LEVEL > 2):
    fp.write('\n pairwise score of this cluster pair is : ' + str(score_val))
    fp.close()
    
  return score_val

##-----------------------------------------------------    
""" this function solves multiple parenting problem (C2)
by uniquely selecting one particular parent
the selection is carried out using a scoring mechanism """
def SolveMultipleParentC2Problem(Output_Text_File):
  for cx in Cluster_Info_Dict:
    
    if (DEBUG_LEVEL > 2):
      fp = open(Output_Text_File, 'a')
      fp.write('\n ***** Examining cluster -- ')
      fp.close()      
      Cluster_Info_Dict[cx]._PrintClusterInfo(cx, Output_Text_File)    
    
    if (Cluster_Info_Dict[cx]._Get_Indegree() > 1):
      # at first form the list to contain the score values for all the child nodes
      # we'll define the score value later
      if (DEBUG_LEVEL > 2):
	fp = open(Output_Text_File, 'a')
	fp.write('\n ***** Examining cluster with more than one indegree -- before in edge list fixing: ')
	fp.close()      
	Cluster_Info_Dict[cx]._PrintClusterInfo(cx, Output_Text_File)
      
      # initialize one dictionary keyed by cluster indices
      scoring_dict = dict()
      for cz in Cluster_Info_Dict[cx]._GetInEdgeList():
	scoring_dict.setdefault(cz, 0)
      # now for each of the in clusters, compute the score 
      for cz in Cluster_Info_Dict[cx]._GetInEdgeList():
	scoring_dict[cz] = ComputeScore(cz, cx, Output_Text_File)
      
      # open the output file
      if (DEBUG_LEVEL > 2):
	fp = open(Output_Text_File, 'a')      
      
      # after computing all such scores for all the in edge clusters
      # we store the values in a list and sort it
      Scoring_List = []
      for cz in Cluster_Info_Dict[cx]._GetInEdgeList():
	if (DEBUG_LEVEL > 2):
	  fp.write('\n scoring dict elem: ' + str(cz) + ' score: ' + str(Scoring_Dict[cz]))
	temp_subl = [cz, scoring_dict[cz]]
	Scoring_List.append(temp_subl)
      
      # after obtaining scores for different taxa set belonging under individual child 
      # nodes of the current node under study, we decide about their ancestral / descendant relationships
      if (DEBUG_LEVEL > 2):
	fp.write('\n --- before sorting the scoring list --- ')
	for i in range(len(Scoring_List)):
	  fp.write('\n elem idx: ' + str(i) + ' cluster label: ' + str(Scoring_List[i][0]) + ' score: ' + str(Scoring_List[i][1]))
	
      # sort the scoring list
      Scoring_List.sort(key=lambda x: x[1])
      if (DEBUG_LEVEL > 2):
	fp.write('\n --- after sorting the scoring list --- ')
	for i in range(len(Scoring_List)):
	  fp.write('\n elem idx: ' + str(i) + ' cluster label: ' + str(Scoring_List[i][0]) + ' score: ' + str(Scoring_List[i][1]))
   
      # close the output file
      if (DEBUG_LEVEL > 2):
	fp.close()
	
      # after sorting the scoring list, now remove all except the last element from the 
      # cx cluster in edge lists
      for i in range(len(Scoring_List) - 1):
	target_delete_clust_idx = Scoring_List[i][0]
	Cluster_Info_Dict[cx]._RemoveInEdge(target_delete_clust_idx)
	Cluster_Info_Dict[target_delete_clust_idx]._RemoveOutEdge(cx)
    
      if (DEBUG_LEVEL > 2):
	fp = open(Output_Text_File, 'a')
	fp.write('\n ***** Examining cluster with more than one indegree -- after in edge list fixing: ')
	fp.close()      
	Cluster_Info_Dict[cx]._PrintClusterInfo(cx, Output_Text_File)
    

##-----------------------------------------------------    
# following code removes extra paranthesis (thus producing insignificant edges) 
# from the tree expression contained in the string
def Remove_Extra_Paranthesis(Final_Supertree_Str):
  L = list(Final_Supertree_Str)
  SL = []	# stack list
  first_bracket_dict = dict()
  
  for i in range(len(L)):
    if (L[i] == '('):
      SL.append(i)
    elif (L[i] == ')'):
      first_bracket_idx = SL.pop()
      first_bracket_dict.setdefault(first_bracket_idx, i)
 
  if (DEBUG_LEVEL > 2):
    print 'L : ', L
    print 'first_bracket_dict: ', first_bracket_dict
  
  for i in range(len(L) - 1):
    if (L[i] == '(') and (L[i+1] == '('):
      sb1 = first_bracket_dict[i]
      sb2 = first_bracket_dict[i+1]
      if (sb1 - sb2 == 1):
	L.pop(i)
	L.insert(i, '$')
	L.pop(sb1)
	L.insert(sb1, '$')
  
  if (DEBUG_LEVEL > 2):
    print 'L : ', L
  
  while (1):
    if '$' in L:
      L.remove('$')
    else:
      break
  
  if (DEBUG_LEVEL > 2):
    print 'L : ', L

  # construct the string containing final supertree
  outstr = ''.join(L)
  return outstr
        
##-----------------------------------------------------        
''' this function returns the root node for the final supertree 
for a depth first forest, multiple root nodes can be possible - so it returns the node with 0 indegree '''
def Extract_Node_Min_Indeg(no_of_clusters):
  min_indeg_node_idx = -1
  valid_node_found = 0
  for i in Cluster_Info_Dict:
    if (Cluster_Info_Dict[i]._GetExploredStatus() == 0):
      if (valid_node_found == 0):
	min_indeg = Cluster_Info_Dict[i]._Get_Indegree()
	min_indeg_node_idx = i
	valid_node_found = 1
      elif (valid_node_found == 1) and (Cluster_Info_Dict[i]._Get_Indegree() < min_indeg):
	min_indeg = Cluster_Info_Dict[i]._Get_Indegree()
	min_indeg_node_idx = i
      elif (valid_node_found == 1) and (Cluster_Info_Dict[i]._Get_Indegree() == min_indeg)\
	    and (Cluster_Info_Dict[i]._Get_Outdegree() > Cluster_Info_Dict[min_indeg_node_idx]._Get_Outdegree()):    
	min_indeg = Cluster_Info_Dict[i]._Get_Indegree()
	min_indeg_node_idx = i
    
  return min_indeg_node_idx
  
##-----------------------------------------------------
# this function prints the tree in Newick format
# sourya - this is the old function with cluster based species list printing in the newick format
def PrintNewick(root_clust_node_idx):
  if 0:
    print 'in function printnewick:   root_clust_node_idx: ', root_clust_node_idx
    print 'taxa set: ', Cluster_Info_Dict[root_clust_node_idx]._GetSpeciesList()  
    print 'out clust list: ', Cluster_Info_Dict[root_clust_node_idx]._GetOutEdgeList()
  
  Tree_Str_List = ''
  # process the node provided it has not been explored yet
  if (Cluster_Info_Dict[root_clust_node_idx]._GetExploredStatus() == 0):  
    # set the explored status of the current node to true
    Cluster_Info_Dict[root_clust_node_idx]._SetExploredStatus()
    # get the out edge list of the current node which are not explored yet 
    outnodes = []
    for l in Cluster_Info_Dict[root_clust_node_idx]._GetOutEdgeList():
      if (Cluster_Info_Dict[l]._GetExploredStatus() == 0):
	outnodes.append(l)
    if (len(outnodes) == 0):
      spec_list = Cluster_Info_Dict[root_clust_node_idx]._GetSpeciesList()
      if (len(spec_list) > 1):
	Tree_Str_List = Tree_Str_List + '('
      Tree_Str_List = Tree_Str_List + ','.join("'" + item + "'" for item in spec_list)
      if (len(spec_list) > 1):
	Tree_Str_List = Tree_Str_List + ')'
    else:
      Tree_Str_List = Tree_Str_List + '('
      Tree_Str_List = Tree_Str_List + ','.join("'" + item + "'" for item in Cluster_Info_Dict[root_clust_node_idx]._GetSpeciesList())
      Tree_Str_List = Tree_Str_List + ','    
      Tree_Str_List = Tree_Str_List + '('
      for i in range(len(outnodes)):
	if (Cluster_Info_Dict[outnodes[i]]._GetExploredStatus() == 0):  
	  Tree_Str_List = Tree_Str_List + PrintNewick(outnodes[i])
	  if (i < (len(outnodes) - 1)):
	    # we check whether any subsequent node belonging to the outnodes list
	    # is left for traverse
	    j = i + 1
	    while (j < len(outnodes)):
	      if (Cluster_Info_Dict[outnodes[j]]._GetExploredStatus() == 0):  
		break
	      j = j + 1	      
	    # in this case, we append one comma
	    if (j < len(outnodes)):
	      Tree_Str_List = Tree_Str_List + ','      
      
      Tree_Str_List = Tree_Str_List + ')'
      Tree_Str_List = Tree_Str_List + ')'
    
  return Tree_Str_List    
    
##-----------------------------------------------------  
''' this function performs transitive reduction of a graph (transitive closure) and subsequently modifies the cluster of nodes
in terms of the edge connectivity, to make it free of redunant edges '''
def CompressDirectedGraph(Reachability_Graph_Mat):    
  no_of_clusters = len(CURRENT_CLUST_IDX_LIST)
  # transitive reduction
  for j in range(no_of_clusters):
    for i in range(no_of_clusters):
      if (Reachability_Graph_Mat[i][j] == 1):
	for k in range(no_of_clusters):
	  if (Reachability_Graph_Mat[j][k] == 1) and (Reachability_Graph_Mat[i][k] == 1):
	    Reachability_Graph_Mat[i][k] = 0
	    # remove the edge from the cluster node directory
	    clust_i = CURRENT_CLUST_IDX_LIST[i]
	    clust_k = CURRENT_CLUST_IDX_LIST[k]
	    Cluster_Info_Dict[clust_i]._RemoveOutEdge(clust_k)
	    Cluster_Info_Dict[clust_k]._RemoveInEdge(clust_i)
    
##-----------------------------------------------------
""" this function creates one new cluster with the given index value
also, it inserts one specified taxa in that cluster """
def Create_Cluster_Taxa_Label(target_clust_idx, target_taxa_label):
  # create the cluster
  Cluster_Info_Dict.setdefault(target_clust_idx, Cluster_node(target_taxa_label))
  # include the cluster idx in the global list CURRENT_CLUST_IDX_LIST
  CURRENT_CLUST_IDX_LIST.append(target_clust_idx)
  # mention the cluster index in the taxa information
  Taxa_Info_Dict[target_taxa_label]._Set_Clust_Idx_taxa_Part(target_clust_idx)
  
##-----------------------------------------------------
""" this function appends one specified taxon on a given cluster """
def Append_Cluster_Taxa_Label(target_clust_idx, target_taxa_label):
  if target_taxa_label not in Cluster_Info_Dict[target_clust_idx]._GetSpeciesList():
    Cluster_Info_Dict[target_clust_idx]._Append_taxa(target_taxa_label)
    # mention the cluster index in the taxa information
    Taxa_Info_Dict[target_taxa_label]._Set_Clust_Idx_taxa_Part(target_clust_idx)  
    
##-----------------------------------------------------
''' this function reads the input tree collection file
the file contains a collection of input candidate source trees
each such tree is composed of a large no of taxa (placed at the leaves of the tree) '''
def Read_Gene_Data_Treelist(ROOTED_TREE, PRESERVE_UNDERSCORE, INPUT_FILE_FORMAT, INPUT_FILENAME):
  ''' depending on the value of INPUT_FILE_FORMAT
  the data is read from the file according to NEWICK or NEXUS format '''
  if (INPUT_FILE_FORMAT == 1):
    Species_TreeList = dendropy.TreeList.get_from_path(INPUT_FILENAME, \
						      schema="newick", \
						      preserve_underscores=PRESERVE_UNDERSCORE, \
						      default_as_rooted=ROOTED_TREE)
  else:
    Species_TreeList = dendropy.TreeList.get_from_path(INPUT_FILENAME, \
						      schema="nexus", \
						      preserve_underscores=PRESERVE_UNDERSCORE, \
						      default_as_rooted=ROOTED_TREE)
  
  return Species_TreeList
    
##-----------------------------------------------------
# this function computes the number of levels between the parent node of node1 and 
# the parent node of node2
# here parent node of node2 is a descendant of parent node of node1
def ComputeLevels(node1, node2):
  lev = 0
  n = node2.parent_node
  while (1):
    if (n == node1.parent_node):
      break
    lev = lev + 1
    if (n.parent_node is not None):
      n = n.parent_node
    else:
      break
  return lev
    
##-----------------------------------------------------
# when node1 and node2 are related via no edge relationship, we compute the total no of levels between them
# this is computed as the sum of distances of these two node levels with their MRCA node
def ComputeLevelsNoEdgeReln(node1, node2, mrca_node, no_parent_prob):
  lev1 = 0
  lev2 = 0
  n1 = node1.parent_node
  n2 = node2.parent_node
  while (1):
    if (n1 == mrca_node):
      break
    lev1 = lev1 + 1
    if (n1.parent_node is not None):
      n1 = n1.parent_node
    else:
      break
    
  while (1):
    if (n2 == mrca_node):
      break
    lev2 = lev2 + 1
    if (n2.parent_node is not None):
      n2 = n2.parent_node
    else:
      break
  
  if (no_parent_prob == 1):
    return (lev1 + lev2)
  else:
    return min(lev1, lev2)

##-----------------------------------------------------
# this function defines relationship between a pair of taxa, with respect to a particular tree
# a taxa is represented by a leaf node in a tree 
# the relationship is either ancestor / descendant, or siblings, or no relationship 
def DefineLeafPairReln(label_taxa1, label_taxa2, Curr_tree, level_info_consider, Cost_update_latest):
  # find the nodes in the tree corresponding to the input taxa
  node1 = Curr_tree.find_node_with_taxon_label(label_taxa1)
  node2 = Curr_tree.find_node_with_taxon_label(label_taxa2)
  # key helps to find out the content in the structure "TaxaPair_Reln_Dict"
  key = (label_taxa1, label_taxa2)
  
  # find the MRCA of these two nodes
  Taxa_Label_List = [label_taxa1, label_taxa2]
  mrca_node = Curr_tree.mrca(taxon_labels=Taxa_Label_List)
      
  if (mrca_node == node1.parent_node) and (mrca_node == node2.parent_node):  #equality relationship 
    TaxaPair_Reln_Dict[key]._AddEdgeCount(BI_DIRECTED_EDGE)
    TaxaPair_Reln_Dict[key]._AddSupportingTree()
    TaxaPair_Reln_Dict[key]._AddLevels(2)
    if (Cost_update_latest == True):
      Taxa_Info_Dict[label_taxa1]._AddOrigEdge(label_taxa2, BI_DIRECTED_EDGE)
      Taxa_Info_Dict[label_taxa2]._AddOrigEdge(label_taxa1, BI_DIRECTED_EDGE)
    if (DEBUG_LEVEL > 2):
      print label_taxa1, ' and ', label_taxa2, 'are connected via BI_DIRECTED_EDGE '
  elif (mrca_node == node1.parent_node):  # checking whether node1 is ancestor of node2
    no_of_levels = ComputeLevels(node1, node2)
    TaxaPair_Reln_Dict[key]._AddSupportingTree()
    TaxaPair_Reln_Dict[key]._AddLevels(no_of_levels)    
    if (level_info_consider == False):    
      TaxaPair_Reln_Dict[key]._AddEdgeCount(DIRECTED_OUT_EDGE)
      if (Cost_update_latest == True):
	Taxa_Info_Dict[label_taxa1]._AddOrigEdge(label_taxa2, DIRECTED_OUT_EDGE)
	Taxa_Info_Dict[label_taxa2]._AddOrigEdge(label_taxa1, DIRECTED_IN_EDGE)
      if (DEBUG_LEVEL > 2):
	print label_taxa1, ' to ', label_taxa2, ' --- DIRECTED_OUT_EDGE '
    else:
      TaxaPair_Reln_Dict[key]._AddEdgeCount(DIRECTED_OUT_EDGE, no_of_levels)
      Taxa_Info_Dict[label_taxa1]._AddOrigEdge(label_taxa2, DIRECTED_OUT_EDGE)
      Taxa_Info_Dict[label_taxa2]._AddOrigEdge(label_taxa1, DIRECTED_IN_EDGE)
      if (DEBUG_LEVEL > 2):
	print label_taxa1, ' to ', label_taxa2, ' --- DIRECTED_OUT_EDGE ', ' no of levels : ', no_of_levels
  elif (mrca_node == node2.parent_node):	# checking whether node2 is ancestor of node1
    no_of_levels = ComputeLevels(node2, node1)
    TaxaPair_Reln_Dict[key]._AddSupportingTree()
    TaxaPair_Reln_Dict[key]._AddLevels(no_of_levels)        
    if (level_info_consider == False):    
      TaxaPair_Reln_Dict[key]._AddEdgeCount(DIRECTED_IN_EDGE)
      if (Cost_update_latest == True):
	Taxa_Info_Dict[label_taxa1]._AddOrigEdge(label_taxa2, DIRECTED_IN_EDGE)
	Taxa_Info_Dict[label_taxa2]._AddOrigEdge(label_taxa1, DIRECTED_OUT_EDGE)
      if (DEBUG_LEVEL > 2):
	print label_taxa1, ' to ', label_taxa2, ' --- DIRECTED_IN_EDGE '
    else:
      TaxaPair_Reln_Dict[key]._AddEdgeCount(DIRECTED_IN_EDGE, no_of_levels)
      if (Cost_update_latest == True):
	Taxa_Info_Dict[label_taxa1]._AddOrigEdge(label_taxa2, DIRECTED_IN_EDGE)
	Taxa_Info_Dict[label_taxa2]._AddOrigEdge(label_taxa1, DIRECTED_OUT_EDGE)
      if (DEBUG_LEVEL > 2):
	print label_taxa1, ' to ', label_taxa2, ' --- DIRECTED_IN_EDGE ', ' no of levels : ', no_of_levels
  else:
    #otherwise the nodes do not have any kind of relationship
    TaxaPair_Reln_Dict[key]._AddSupportingTree()
    no_of_levels = ComputeLevelsNoEdgeReln(node1, node2, mrca_node, 1)
    TaxaPair_Reln_Dict[key]._AddLevels(no_of_levels)        
    if (level_info_consider == False):    
      TaxaPair_Reln_Dict[key]._AddEdgeCount(NO_EDGE)
      if (Cost_update_latest == True):
	Taxa_Info_Dict[label_taxa1]._AddOrigEdge(label_taxa2, NO_EDGE)
	Taxa_Info_Dict[label_taxa2]._AddOrigEdge(label_taxa1, NO_EDGE)  
      if (DEBUG_LEVEL > 2):
	print label_taxa1, ' to ', label_taxa2, ' --- NO_EDGE '    
    else:
      no_of_levels = ComputeLevelsNoEdgeReln(node1, node2, mrca_node, 0)
      TaxaPair_Reln_Dict[key]._AddEdgeCount(NO_EDGE, no_of_levels)
      if (Cost_update_latest == True):
	Taxa_Info_Dict[label_taxa1]._AddOrigEdge(label_taxa2, NO_EDGE)
	Taxa_Info_Dict[label_taxa2]._AddOrigEdge(label_taxa1, NO_EDGE)    
      if (DEBUG_LEVEL > 2):
	print label_taxa1, ' to ', label_taxa2, ' --- NO_EDGE ', ' no of levels : ', no_of_levels    


#--------------------------------------------------
# this function returns the label of an internal or a leaf node 
# in terms of newick representation
def Node_Label(inp_node):
  return str(inp_node.as_newick_string(suppress_edge_lengths=True))
