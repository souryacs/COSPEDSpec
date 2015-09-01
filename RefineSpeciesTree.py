#!/usr/bin/env python

import Header
from Header import *
import UtilFunc
from UtilFunc import *

#-------------------------------------------
# this function computes the level score corresponding to the cluster pair provided
# for 1st method of binary refinement
def ComputeLevelScore(clust_spec_list1, clust_spec_list2):
  level_score = 0
  for x in clust_spec_list1:
    for y in clust_spec_list2:
      key1 = (x, y)
      key2 = (y, x)
      if (key1 in TaxaPair_Reln_Dict):
	level_score = level_score + TaxaPair_Reln_Dict[key1]._GetAvgLevelPerTree()
      elif (key2 in TaxaPair_Reln_Dict):
	level_score = level_score + TaxaPair_Reln_Dict[key2]._GetAvgLevelPerTree()

  return level_score

#-------------------------------------------
# this function processes one internal node (basically the children list)
# to resolve multifurcation
def ResolveMultifurcation(Gene_TreeList, Curr_tree, clust_species_list, no_of_input_clusters, Output_Text_File, Bin_Sup_Opt):
  # total number of clusters
  no_of_clust = no_of_input_clusters
  
  # allocate a 2D square matrix of no_of_clust dimension
  Level_Score_Mat = numpy.zeros((no_of_clust, no_of_clust), dtype=numpy.int)
  Norm_Level_Score_Mat = numpy.zeros((no_of_clust, no_of_clust), dtype=numpy.float)
    
  if (Bin_Sup_Opt == 1):      
    # compute the pairwise score of clusters
    # according to the mutual level information
    # level value is computed by comparing pairwise taxa level information 
    # from individual species list
    for i in range(no_of_clust - 1):
      for j in range(i+1, no_of_clust):
    	if (DEBUG_LEVEL > 2):
	  fp = open(Output_Text_File, 'a')
	  fp.write('\n taxa list pair analysis:  i - ' + str(i) + '  j - ' + str(j) + \
	    ' 1st list : ' + str(clust_species_list[i]) + ' 2nd list - ' + str(clust_species_list[j]))	
	  fp.close() 
	# this function computes the score of the cluster pair with level wise information
	level_score = ComputeLevelScore(clust_species_list[i], clust_species_list[j])    
	# insert the score in a matrix
	Level_Score_Mat[i][j] = level_score
	Level_Score_Mat[j][i] = level_score
  else:
    # here we compute the ILS score of the current cluster pair
    # with respect to input gene tree list
    if (DEBUG_LEVEL > 2):
      fp = open(Output_Text_File, 'a')
      fp.write('\n Examining ILS score for individual cluster pairs ')
      fp.close()      
    for tr in Gene_TreeList:
      Inp_Tree = Tree(tr)
      for i in range(no_of_clust - 1):
	for j in range(i+1, no_of_clust):
	  Clust_Pair_Taxa_All = []
	  for x in clust_species_list[i]:
	    Clust_Pair_Taxa_All.append(x)
	  for x in clust_species_list[j]:
	    Clust_Pair_Taxa_All.append(x)		  
	  # if all taxa under analysis is present in the input gene tree
	  # then find the MRCA of the taxa list in the input gene tree
	  MRCA_Node_Taxa_List = Inp_Tree.mrca(taxon_labels=Clust_Pair_Taxa_All)
	  if (MRCA_Node_Taxa_List is not None):
	    XL_MRCA_node = len(MRCA_Node_Taxa_List.leaf_nodes()) - len(Clust_Pair_Taxa_All)
	    Level_Score_Mat[i][j] = Level_Score_Mat[i][j] + XL_MRCA_node
	    Level_Score_Mat[j][i] = Level_Score_Mat[j][i] + XL_MRCA_node
         
  # this variable signifies whether we need to check the leaf nodes 
  # for clustering
  LEAF_NODE_CHECK_REQUIRED = True
        
  # loop to execute the agglomerative clustering
  while(no_of_clust > 2): 
    
    # add - sourya
    flag = 0
    if (LEAF_NODE_CHECK_REQUIRED == True):
      # this check will be done at the very first iteration
      # provided that we follow MRCA based refinement
      if 1:	#(no_of_clust == no_of_input_clusters) and (Bin_Sup_Opt == 2):
	min_idx_i = 0
	min_idx_j = 1      
	for i in range(len(clust_species_list) - 1):
	  if (flag == 1):
	    break
	  if (len(clust_species_list[i]) == 1):
	    for j in range(i+1, len(clust_species_list)):
	      if (len(clust_species_list[j]) == 1):
		# check if the species pair has final selected edge of type 0
		spec1 = clust_species_list[i][0]
		spec2 = clust_species_list[j][0]
		# check whether the taxa pair was selected for bi directed edge relation
		if (Taxa_Info_Dict[spec1]._Get_Taxa_Part_Clust_Idx() == Taxa_Info_Dict[spec2]._Get_Taxa_Part_Clust_Idx()):
		  min_idx_i = i
		  min_idx_j = j
		  flag = 1
		  break
    # end add - sourya
    
    if (DEBUG_LEVEL > 2):
      fp = open(Output_Text_File, 'a')
      fp.write('\n iteration start --- number of clusters: ' + str(no_of_clust))
      fp.write('\n clust_species_list : ' + str(clust_species_list))
      fp.write('\n printing contents of Level_Score_Mat ---- ')
      for i in range(no_of_clust):
	fp.write('\n ')
	for j in range(no_of_clust):
	  fp.write(' ' + str(Level_Score_Mat[i][j]))     
      fp.close()
    
    # allocate one new square matrix which will contain the 
    # normalized matrix elements (w.r.t the sum of sum of rows and columns)
    sum_list = []
    for i in range(no_of_clust):
      t = 0
      for j in range(no_of_clust):
	t = t + Level_Score_Mat[i][j]
      sum_list.append(t)

    for i in range(no_of_clust - 1):
      for j in range(i+1, no_of_clust):
	# comment - sourya
	# previously we normalized the extra lineage based score
	# by the sum of extra lineages for all other taxa from the taxa indices i and j
	Norm_Level_Score_Mat[i][j] = (Level_Score_Mat[i][j] * 1.0) / (sum_list[i] + sum_list[j])
	# add - sourya
	# here we do not normalize the extra lineages although we keep the same naming convention
	# here we check the following quantity to be lowest as possible
	#Norm_Level_Score_Mat[i][j] = (3 * Level_Score_Mat[i][j] - sum_list[i] - sum_list[j])
	# end add - sourya
	Norm_Level_Score_Mat[j][i] = Norm_Level_Score_Mat[i][j]
    
    if (DEBUG_LEVEL > 2):
      fp = open(Output_Text_File, 'a')
      fp.write('\n printing contents of sum_list --- ' + str(sum_list))
      fp.write('\n printing contents of Norm_Level_Score_Mat ---- ')
      for i in range(no_of_clust):
	fp.write('\n ')
	for j in range(no_of_clust):
	  fp.write(' ' + str(Norm_Level_Score_Mat[i][j]))   
      fp.close()
    
    if (flag == 0): #or (Bin_Sup_Opt == 1):	# this condition is added - sourya
      LEAF_NODE_CHECK_REQUIRED = False
      min_val = Norm_Level_Score_Mat[0][1]
      min_idx_i = 0
      min_idx_j = 1
      for i in range(no_of_clust - 1):
	for j in range(i+1, no_of_clust):
	  if (i == j):
	    continue
	  if (Norm_Level_Score_Mat[i][j] < min_val):
	    min_val = Norm_Level_Score_Mat[i][j]
	    min_idx_i = i
	    min_idx_j = j

    # note down the taxa list in these two indices of the clust_species_list
    taxa_list = []
    for x in clust_species_list[min_idx_i]:
      taxa_list.append(x)
    for x in clust_species_list[min_idx_j]:
      taxa_list.append(x)

    if (DEBUG_LEVEL > 2):
      fp = open(Output_Text_File, 'a')
      fp.write('\n min_idx_i ' + str(min_idx_i) + ' min_idx_j : ' + str(min_idx_j))
      fp.write('\n min_idx_i species list ' + str(clust_species_list[min_idx_i]))
      fp.write('\n min_idx_j species list ' + str(clust_species_list[min_idx_j]))
      fp.write('\n complete species list ' + str(taxa_list))
      fp.close()
      
    #---------------------------------------------------------      
    # for individual clusters, we check if the cluster contains one or more species
    # case 1 - both the clusters have > 1 species
    # and the clusters are represented by an internal node which is the MRCA of the constituent species set
    if (len(clust_species_list[min_idx_i]) > 1) and (len(clust_species_list[min_idx_j]) > 1):
      first_cluster_mrca_node = Curr_tree.mrca(taxon_labels=clust_species_list[min_idx_i])
      second_cluster_mrca_node = Curr_tree.mrca(taxon_labels=clust_species_list[min_idx_j])
      all_taxa_mrca_node = Curr_tree.mrca(taxon_labels=taxa_list)
      if (DEBUG_LEVEL > 2):
	fp = open(Output_Text_File, 'a')
	fp.write('\n label of first_cluster_mrca_node: ' + str(Node_Label(first_cluster_mrca_node)))      
	fp.write('\n label of second_cluster_mrca_node: ' + str(Node_Label(second_cluster_mrca_node)))
	fp.write('\n label of all_taxa_mrca_node: ' + str(Node_Label(all_taxa_mrca_node)))
	fp.close()
	
      # create new internal node 
      newnode = Node()
      # its parent node will be the previous MRCA node of all the taxa in two clusters
      all_taxa_mrca_node.add_child(newnode)
      newnode.parent_node = all_taxa_mrca_node
      all_taxa_mrca_node.remove_child(first_cluster_mrca_node)
      first_cluster_mrca_node.parent_node = None
      all_taxa_mrca_node.remove_child(second_cluster_mrca_node)
      second_cluster_mrca_node.parent_node = None
      # add these individual clusters' MRCA node as its children
      newnode.add_child(first_cluster_mrca_node)
      first_cluster_mrca_node.parent_node = newnode
      newnode.add_child(second_cluster_mrca_node)
      second_cluster_mrca_node.parent_node = newnode
      # update splits of the resulting tree
      Curr_tree.update_splits(delete_outdegree_one=True)
      
    # case 2 and 3 - one cluster has at least 2 species, while other is a leaf
    elif (len(clust_species_list[min_idx_i]) == 1) and (len(clust_species_list[min_idx_j]) > 1):
      first_cluster_leaf_node = Curr_tree.find_node_with_taxon_label(clust_species_list[min_idx_i][0])
      second_cluster_mrca_node = Curr_tree.mrca(taxon_labels=clust_species_list[min_idx_j])
      all_taxa_mrca_node = Curr_tree.mrca(taxon_labels=taxa_list)
      if (DEBUG_LEVEL > 2):
	fp = open(Output_Text_File, 'a')
	fp.write('\n first cluster is a leaf - its label: ' + str(Node_Label(first_cluster_leaf_node)))      
	fp.write('\n label of second_cluster_mrca_node: ' + str(Node_Label(second_cluster_mrca_node)))
	fp.write('\n label of all_taxa_mrca_node: ' + str(Node_Label(all_taxa_mrca_node)))
	fp.close()
      
      # create new internal node 
      newnode = Node()
      # its parent node will be the previous MRCA node of all the taxa in two clusters
      all_taxa_mrca_node.add_child(newnode)
      newnode.parent_node = all_taxa_mrca_node
      all_taxa_mrca_node.remove_child(first_cluster_leaf_node)
      first_cluster_leaf_node.parent_node = None
      all_taxa_mrca_node.remove_child(second_cluster_mrca_node)
      second_cluster_mrca_node.parent_node = None
      # add these individual clusters' MRCA node as its children
      newnode.add_child(first_cluster_leaf_node)
      first_cluster_leaf_node.parent_node = newnode
      newnode.add_child(second_cluster_mrca_node)
      second_cluster_mrca_node.parent_node = newnode
      # update splits of the resulting tree
      Curr_tree.update_splits(delete_outdegree_one=True)
      
    elif (len(clust_species_list[min_idx_i]) > 1) and (len(clust_species_list[min_idx_j]) == 1):
      first_cluster_mrca_node = Curr_tree.mrca(taxon_labels=clust_species_list[min_idx_i])
      second_cluster_leaf_node = Curr_tree.find_node_with_taxon_label(clust_species_list[min_idx_j][0])
      all_taxa_mrca_node = Curr_tree.mrca(taxon_labels=taxa_list)
      if (DEBUG_LEVEL > 2):
	fp = open(Output_Text_File, 'a')
	fp.write('\n label of first_cluster_mrca_node: ' + str(Node_Label(first_cluster_mrca_node)))      
	fp.write('\n label of second_cluster_mrca_node: ' + str(Node_Label(second_cluster_leaf_node)))
	fp.write('\n label of all_taxa_mrca_node: ' + str(Node_Label(all_taxa_mrca_node)))
	fp.close()
      
      # create new internal node 
      newnode = Node()
      # its parent node will be the previous MRCA node of all the taxa in two clusters
      all_taxa_mrca_node.add_child(newnode)
      newnode.parent_node = all_taxa_mrca_node
      all_taxa_mrca_node.remove_child(first_cluster_mrca_node)
      first_cluster_mrca_node.parent_node = None
      all_taxa_mrca_node.remove_child(second_cluster_leaf_node)
      second_cluster_leaf_node.parent_node = None
      # add these individual clusters' MRCA node as its children
      newnode.add_child(first_cluster_mrca_node)
      first_cluster_mrca_node.parent_node = newnode
      newnode.add_child(second_cluster_leaf_node)
      second_cluster_leaf_node.parent_node = newnode
      # update splits of the resulting tree
      Curr_tree.update_splits(delete_outdegree_one=True)
      
    # case 4 - when both child clusters are leaf nodes 
    else:
      first_cluster_leaf_node = Curr_tree.find_node_with_taxon_label(clust_species_list[min_idx_i][0])
      second_cluster_leaf_node = Curr_tree.find_node_with_taxon_label(clust_species_list[min_idx_j][0])
      all_taxa_mrca_node = Curr_tree.mrca(taxon_labels=taxa_list)
      if (DEBUG_LEVEL > 2):
	fp = open(Output_Text_File, 'a')
	fp.write('\n first cluster is a leaf - its label: ' + str(Node_Label(first_cluster_leaf_node)))      
	fp.write('\n second cluster is a leaf - its label: ' + str(Node_Label(second_cluster_leaf_node)))
	fp.write('\n label of all_taxa_mrca_node: ' + str(Node_Label(all_taxa_mrca_node)))
	fp.close()
      
      # create new internal node 
      newnode = Node()
      # its parent node will be the previous MRCA node of all the taxa in two clusters
      all_taxa_mrca_node.add_child(newnode)
      newnode.parent_node = all_taxa_mrca_node
      all_taxa_mrca_node.remove_child(first_cluster_leaf_node)
      first_cluster_leaf_node.parent_node = None
      all_taxa_mrca_node.remove_child(second_cluster_leaf_node)
      second_cluster_leaf_node.parent_node = None
      # add these individual clusters' MRCA node as its children
      newnode.add_child(first_cluster_leaf_node)
      first_cluster_leaf_node.parent_node = newnode
      newnode.add_child(second_cluster_leaf_node)
      second_cluster_leaf_node.parent_node = newnode
      # update splits of the resulting tree
      Curr_tree.update_splits(delete_outdegree_one=True)
    
    #---------------------------------------------------------      
    if (DEBUG_LEVEL > 2):
      fp = open(Output_Text_File, 'a')
      fp.write('\n label of newnode: ' + str(Node_Label(newnode)))
      fp.write('\n label of all taxa mrca node (recomputed): ' + str(Node_Label(Curr_tree.mrca(taxon_labels=taxa_list))))      
      fp.write('\n before inserting row col, Level_Score_Mat dimension: ' + str(Level_Score_Mat.size))
      fp.close()
    
    # adjust the Level_Score_Mat by inserting one new row and column corresponding to the new cluster
    # and then deleting the information of earlier two clusters
    # first append one row
    Level_Score_Mat = numpy.vstack((Level_Score_Mat, numpy.zeros((1, no_of_clust), dtype=numpy.int)))
    # then append one column
    Level_Score_Mat = numpy.hstack((Level_Score_Mat, numpy.zeros((no_of_clust + 1, 1), dtype=numpy.int)))
    
    if (DEBUG_LEVEL > 2):
      fp = open(Output_Text_File, 'a')
      fp.write('\n after inserting row col, Level_Score_Mat dimension: ' + str(Level_Score_Mat.size))
      fp.close()
      
    Level_Score_Mat = numpy.reshape(Level_Score_Mat, ((no_of_clust + 1), (no_of_clust + 1)), order='C')
    
    if (Bin_Sup_Opt == 1):
      # now fill the elements of the new added row and column
      # it is the average of elements located corresponding to min_idx_i and min_idx_j
      for k in range(no_of_clust):
	Level_Score_Mat[k][no_of_clust] = (Level_Score_Mat[k][min_idx_i] + Level_Score_Mat[k][min_idx_j]) / 2 - 1
	Level_Score_Mat[no_of_clust][k] = Level_Score_Mat[k][no_of_clust]
    else:
      # MRCA based adjustment
      # now fill the elements of the new added row and column
      for k in range(no_of_clust):
	# for any index k, the number of extra lineage is the maximum of these three quantities
	Level_Score_Mat[k][no_of_clust] = max(Level_Score_Mat[k][min_idx_i], Level_Score_Mat[k][min_idx_j], \
	  Level_Score_Mat[min_idx_i][min_idx_j])
	Level_Score_Mat[no_of_clust][k] = Level_Score_Mat[k][no_of_clust]
    
    # now remove the rows and columns corresponding to min_idx_i and min_idx_j
    Level_Score_Mat = numpy.delete(Level_Score_Mat, (min_idx_i), axis=0)	# delete the row
    Level_Score_Mat = numpy.delete(Level_Score_Mat, (min_idx_i), axis=1)	# delete the column
    Level_Score_Mat = numpy.delete(Level_Score_Mat, (min_idx_j - 1), axis=0)	# delete the row
    Level_Score_Mat = numpy.delete(Level_Score_Mat, (min_idx_j - 1), axis=1)	# delete the column

    # clear Norm_Level_Score_Mat
    Norm_Level_Score_Mat = numpy.delete(Norm_Level_Score_Mat, (min_idx_i), axis=0)	# delete the row
    Norm_Level_Score_Mat = numpy.delete(Norm_Level_Score_Mat, (min_idx_i), axis=1)	# delete the column
    Norm_Level_Score_Mat.fill(0)
    
    # remove individual clusters' taxa information from the clust_species_list
    # and add taxa_list as a new element
    clust_species_list.pop(min_idx_i)
    clust_species_list.pop(min_idx_j - 1)
    clust_species_list.append(taxa_list)    
    
    # decrement the number of clusters considered
    no_of_clust = no_of_clust - 1
        
    #if (Bin_Sup_Opt == 2) and (no_of_clust > 2):    
      ## fill the elements of the last row and column of the Level_Score_Mat
      #j = no_of_clust - 1
      #for i in range(j):
	#if (DEBUG_LEVEL > 2):
	  #fp = open(Output_Text_File, 'a')
	  #fp.write('\n Examining ILS score for the clusters ' + str(i) + ' and ' + str(j) + \
	    #' species list ' + str(clust_species_list[i]) + ' and ' + str(clust_species_list[j]))
	  #fp.close()
	
	## this function computes the ILS score of the current cluster pair
	## with respect to input gene tree list
	#level_score = ComputeDCScore(Gene_TreeList, clust_species_list[i], clust_species_list[j], Output_Text_File)
	#Level_Score_Mat[i][j] = level_score
	#Level_Score_Mat[j][i] = level_score
      
  return
            
#-------------------------------------------
# this function refines input supertree such that the supertree becomes binary
# this is required for proper benchmarking with existing binary tree construction methods on 
# ILS sorting
def Refine_Supertree_Binary_Form(Gene_TreeList, Curr_tree, Output_Text_File, Bin_Sup_Opt):
  # we traverse input tree internal nodes in postorder fashion
  # and list the child nodes of it
  # if the no of children > 2 then it is a case of multifurcation
  # for resolving
  for curr_node in Curr_tree.postorder_internal_node_iter():
    curr_node_children = curr_node.child_nodes()
    if (len(curr_node_children) > 2):
      # create a list which will contain the species list lying under 
      # individual child nodes of rhe current node
      clust_species_list = []
      # examine individual nodes of the current node's children list
      for x in curr_node_children:
	if (x.is_leaf() == True):
	  subl = []
	  subl.append(x.taxon.label)
	else:
	  subl = GetTaxaUnderInternalNode(x)
	clust_species_list.append(subl)
      
      # call the resolving routine
      ResolveMultifurcation(Gene_TreeList, Curr_tree, clust_species_list, len(curr_node_children), Output_Text_File, Bin_Sup_Opt)

#-----------------------------------------------------
# this is the taxa list generated from current internal node
def GetTaxaUnderInternalNode(curr_node):
  taxa_list_from_curr_internal_node = []
  for n in curr_node.leaf_nodes():
    taxa_list_from_curr_internal_node.append(n.taxon.label)
  return taxa_list_from_curr_internal_node


