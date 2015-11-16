#!/usr/bin/env python

import Header
from Header import *
import UtilFunc
from UtilFunc import *

#-------------------------------------------
"""
this function fills the distance matrix using normalized excess gene count
for a particular taxa cluster, it uses one representative taxon of that taxa cluster
"""
def Fill_DistMat_SingleEntry(DistMat, no_of_clust, clust_species_list):
	"""
	check and explore each pair of taxa clusters
	"""
	for i in range(no_of_clust - 1):
		for j in range(i+1, no_of_clust):
			# here both clust_species_list[i] and clust_species_list[j]
			# are one element lists (according to their construction)
			# we have extracted the corresponding element by using [0] operator (extracting first element)
			elem_found = False
			for k1 in range(len(clust_species_list[i])):
				for k2 in range(len(clust_species_list[j])):  
					x1 = clust_species_list[i][k1]
					x2 = clust_species_list[j][k2]
					key1 = (x1, x2)
					key2 = (x2, x1)
					#print 'key1: ', key1, ' key2: ', key2
					if key1 in TaxaPair_Reln_Dict:
						elem_found = True
						DistMat[j][i] = DistMat[i][j] = TaxaPair_Reln_Dict[key1]._GetNormalizedXLSumGeneTrees()
						break
					elif key2 in TaxaPair_Reln_Dict:
						elem_found = True
						DistMat[j][i] = DistMat[i][j] = TaxaPair_Reln_Dict[key2]._GetNormalizedXLSumGeneTrees()
						break
				if (elem_found == True):
					break
			if (elem_found == False):
				DistMat[j][i] = DistMat[i][j] = 0

	return
	
#-------------------------------------------
"""
this function fills the distance matrix using normalized excess gene count
for a particular taxa cluster, it uses aerage information of that taxa cluster
"""
def Fill_DistMat_AvgEntry(DistMat, no_of_clust, clust_species_list):
	"""
	check and explore each pair of taxa clusters
	"""
	for i in range(no_of_clust - 1):
		for j in range(i+1, no_of_clust):
			# here both clust_species_list[i] and clust_species_list[j]
			# are one element lists (according to their construction)
			# we have extracted the corresponding element by using [0] operator (extracting first element)
			curr_taxa_pair_list = []
			for k1 in range(len(clust_species_list[i])):
				for k2 in range(len(clust_species_list[j])):	  
					x1 = clust_species_list[i][k1]
					x2 = clust_species_list[j][k2]
					key1 = (x1, x2)
					key2 = (x2, x1)
					#print 'key1: ', key1, ' key2: ', key2
					if key1 in TaxaPair_Reln_Dict:
						curr_taxa_pair_list.append(TaxaPair_Reln_Dict[key1]._GetNormalizedXLSumGeneTrees())
					else:
						curr_taxa_pair_list.append(TaxaPair_Reln_Dict[key2]._GetNormalizedXLSumGeneTrees())
			# average of this pairwise list is used as the XT approximation
			if (len(curr_taxa_pair_list) > 0):
				DistMat[j][i] = DistMat[i][j] = (sum(curr_taxa_pair_list) * 1.0) / len(curr_taxa_pair_list)
			else:
				DistMat[j][i] = DistMat[i][j] = 0
	
	return

#-------------------------------------------
"""
this function finds a single minimum from the input matrix
"""
def Find_Unique_Min(Norm_DistMat, no_of_clust, clust_species_list):
	# initialize
	min_val = Norm_DistMat[0][1]
	min_idx_i = 0
	min_idx_j = 1
	# traverse through the matrix elements
	for i in range(no_of_clust - 1):
		for j in range(i+1, no_of_clust):
			if (i == j):
				continue
			if (Norm_DistMat[i][j] < min_val):
				min_val = Norm_DistMat[i][j]
				min_idx_i = i
				min_idx_j = j
			elif (Norm_DistMat[i][j] == min_val):
				# comment - sourya
				## here we prioritize the cluster pair having minimum number of species
				#if (len(clust_species_list[i]) + len(clust_species_list[j])) < (len(clust_species_list[min_idx_i]) + len(clust_species_list[min_idx_j])):
					#min_idx_i = i
					#min_idx_j = j
				# end comment - sourya
				# add - sourya
				if (Get_R1R2Reln_ValDiff_Abs(clust_species_list[i][0], clust_species_list[j][0]) < \
					Get_R1R2Reln_ValDiff_Abs(clust_species_list[min_idx_i][0], clust_species_list[min_idx_j][0])):
					min_idx_i = i
					min_idx_j = j
				# end add - sourya
	
	return min_idx_i, min_idx_j

#--------------------------------------------------------
# this function is a shortcut to obtain the normalized expression 
# used in the agglomerative clustering proposed in this code
# as various methods are experimented, corresponding various forms of 
# agglomerative clustering is tried
#--------------------------------------------------------
def ObtainNormalizedVal(num, denom1, denom2):
  if ((denom1 + denom2) > 0):
    return (num * 1.0) / (denom1 + denom2)
  else:
    return 0

##---------------------------------------------
""" 
function to print the matrix content
N is the matrix dimension
"""
def PrintMatrixContent(N, TaxaList, inp_data, inp_str, textfile):
  fp = open(textfile, 'a')
  fp.write('\n printing contents of ' + str(inp_str) + ' ---- ')
  for i in range(N):
    fp.write('\n ' + str(i) + '--' + str(TaxaList[i]) + '--->>')
    for j in range(i+1):
      fp.write(' ' + str(inp_data[i][j]))
  fp.close()


#-------------------------------------------------------------
"""
this function checks whether the input pair of taxa list is basically 
a couplet, belonging to an existing set of siblings
"""
def CheckEstablishedSibling(taxa_list1, taxa_list2):
	if (len(taxa_list1) == 1) and (len(taxa_list2) == 1):
		coup1 = [taxa_list1[0], taxa_list2[0]]
		coup2 = [taxa_list2[0], taxa_list1[0]]
		#print 'coup1: ' + str(coup1) + '  or coup2: ' + str(coup2) + ' is an established sibling pair'
		if (coup1 in Sibling_Couplet_List) or (coup2 in Sibling_Couplet_List):
			return True
	return False
	
#-----------------------------------------------------------------
def Get_R1R2Reln_ValDiff(x1, x2):
	key1 = (x1, x2)
	key2 = (x2, x1)
	if key1 in TaxaPair_Reln_Dict:	
		return TaxaPair_Reln_Dict[key1]._GetR1R2LevelDiff()
	elif key2 in TaxaPair_Reln_Dict:
		return (-1) * TaxaPair_Reln_Dict[key2]._GetR1R2LevelDiff()
	
	return 0

#-----------------------------------------------------------------
def Get_R1R2Reln_ValDiff_Abs(x1, x2):
	key1 = (x1, x2)
	key2 = (x2, x1)
	if key1 in TaxaPair_Reln_Dict:	
		return TaxaPair_Reln_Dict[key1]._GetR1R2AbsLevelDiff()
	elif key2 in TaxaPair_Reln_Dict:
		return TaxaPair_Reln_Dict[key2]._GetR1R2AbsLevelDiff()
	
	return 0

#-----------------------------------------------------------------
"""
this function checks a couplet whether they can be related by R1 relation
according to the all relation level count difference statistic
"""
def Check_R1Reln_Majority(x1, x2):
	key1 = (x1, x2)
	key2 = (x2, x1)
	if key1 in TaxaPair_Reln_Dict:
		if (TaxaPair_Reln_Dict[key1]._CheckR2RelnLevelConsensus()):
			return -1
		elif (TaxaPair_Reln_Dict[key1]._CheckR1RelnLevelConsensus()):
			return 1
	elif key2 in TaxaPair_Reln_Dict:
		if (TaxaPair_Reln_Dict[key2]._CheckR2RelnLevelConsensus()):
			return 1
		elif (TaxaPair_Reln_Dict[key2]._CheckR1RelnLevelConsensus()):
			return -1
	
	return KEY_ABSENCE_INDICATOR	# key absence indicator
	
#-----------------------------------------------------------------
"""
this function checks according to the all relation level count difference statistic
which cluster can be ancestor and which can be descendant
returns - 1) binary value 1 or 0 depending on the couplet is a sibling or not
2) indices x and y where the species in index x is placed at a higher level than the species at index y
with respect to the input gene trees' configurations
"""
def CheckR1RelationLevelBased(x1, x2):
	key1 = (x1, x2)
	key2 = (x2, x1)
	if key1 in TaxaPair_Reln_Dict:
		if (TaxaPair_Reln_Dict[key1]._CheckR3RelnLevelConsensus()):
			return 0
		elif (TaxaPair_Reln_Dict[key1]._CheckHigherR2RelnLevelValue()):
			return -1
		elif (TaxaPair_Reln_Dict[key1]._CheckHigherR1RelnLevelValue()):
			return 1
	elif key2 in TaxaPair_Reln_Dict:
		if (TaxaPair_Reln_Dict[key2]._CheckR3RelnLevelConsensus()):
			return 0
		elif (TaxaPair_Reln_Dict[key2]._CheckHigherR2RelnLevelValue()):
			return 1
		elif (TaxaPair_Reln_Dict[key2]._CheckHigherR1RelnLevelValue()):
			return -1
	
	return KEY_ABSENCE_INDICATOR	# key absence indicator

#-------------------------------------------
"""
this function processes input distance matrix in every iteration
and finds the pair of indices satisfying minimum distance criterion 
used in NJ based algorithm
"""
def Get_NJ_Based_Min_Pair_Idx(DistMat, Norm_DistMat, no_of_clust, clust_species_list, NJ_RULE_USED, Output_Text_File):
	
	if (DEBUG_LEVEL >= 2):
		fp = open(Output_Text_File, 'a')
		fp.write('\n\n\n\n\n\n **** iteration start --- number of clusters: ' + str(no_of_clust))
		fp.write('\n clust_species_list : ' + str(clust_species_list))
		fp.close()
		PrintMatrixContent(no_of_clust, clust_species_list, DistMat, 'DistMat', Output_Text_File)
	
	# allocate one new square matrix which will contain the 
	# normalized matrix elements (w.r.t the sum of sum of rows and columns)
	sum_list = []
	for i in range(no_of_clust):
		t = 0
		for j in range(no_of_clust):
			t = t + DistMat[i][j]
		sum_list.append(t)

	for i in range(no_of_clust - 1):
		for j in range(i+1, no_of_clust):
			if (NJ_RULE_USED == AGGLO_CLUST):
				# we normalize the extra lineage based score
				# by the sum of extra lineages for all other taxa from the taxa indices i and j
				# modified - sourya
				#Norm_DistMat[i][j] = ObtainNormalizedVal(DistMat[i][j], sum_list[i], sum_list[j])
				# add - sourya
				Norm_DistMat[i][j] = ObtainNormalizedVal(DistMat[i][j], sum_list[i] - DistMat[i][j], sum_list[j] - DistMat[i][j])
				# end add - sourya
				Norm_DistMat[j][i] = Norm_DistMat[i][j]
			else:
				ri = sum_list[i] / (no_of_clust - 2)
				rj = sum_list[j] / (no_of_clust - 2)
				Norm_DistMat[i][j] = (DistMat[i][j] - ri - rj)
				Norm_DistMat[j][i] = Norm_DistMat[i][j]
				
	if (DEBUG_LEVEL >= 2):
		fp = open(Output_Text_File, 'a')
		fp.write('\n printing contents of sum_list --- ' + str(sum_list))
		fp.close()
		PrintMatrixContent(no_of_clust, clust_species_list, Norm_DistMat, 'Norm_DistMat', Output_Text_File)

	# now we have to find the minimum among these elements 
	# present in the matrix Norm_DistMat
	min_idx_i, min_idx_j = Find_Unique_Min(Norm_DistMat, no_of_clust, clust_species_list)

	return min_idx_i, min_idx_j

#-------------------------------------------
"""
checks whether a taxa cluster specified by the input index is a leaf
"""
def IsLeafCluster(clust_species_list, idx):
	if (len(clust_species_list[idx]) == 1):
		return True
	return False

#-------------------------------------------
"""
this function places one subtree at an edge of a second subtree
also adjusts their parent information
parameters:
1) source subtree (which needs to be re-positioned)
2) edge of destination subtree - indicated by the child node
3) Curr_tree: Tree containing all these subtrees

It creates one new internal node within the edge of the destination subtree
and places the source subtree as its children
"""
def InsertSubTree(Curr_tree, src_subtree, child_dest_subtree, Output_Text_File):
	
	if (DEBUG_LEVEL >= 2):
		fp = open(Output_Text_File, 'a')
		fp.write('\n label of src_subtree: ' + str(Node_Label(src_subtree)))      
		fp.write('\n label of child_dest_subtree: ' + str(Node_Label(child_dest_subtree)))
		fp.close()
	
	# create new internal node 
	newnode = Node()
	
	parent_dest_subtree = child_dest_subtree.parent_node
	parent_src_subtree = src_subtree.parent_node
	
	# its parent node will be the previous MRCA node of all the taxa in two clusters
	parent_dest_subtree.add_child(newnode)
	newnode.parent_node = parent_dest_subtree
	parent_dest_subtree.remove_child(child_dest_subtree)
	child_dest_subtree.parent_node = None
	newnode.add_child(child_dest_subtree)
	child_dest_subtree.parent_node = newnode
	
	if (parent_src_subtree is not None):
		parent_src_subtree.remove_child(src_subtree)
		src_subtree.parent_node = None
	newnode.add_child(src_subtree)
	src_subtree.parent_node = newnode
	
	# update splits of the resulting tree
	Curr_tree.update_splits(delete_outdegree_one=False)

	if (DEBUG_LEVEL >= 2):
		fp = open(Output_Text_File, 'a')
		fp.write('\n label of newnode: ' + str(Node_Label(newnode)))
		fp.write('\n label of dest_subtree: ' + str(Node_Label(child_dest_subtree.parent_node)))  
		fp.close()

	return Curr_tree

#-------------------------------------------
"""
this function has following parameters:
1) first_cluster_mrca_node: root of 1st subtree 
2) second_cluster_mrca_node: root of 2nd subtree 
3) all_taxa_mrca_node: root of all these trees
4) Curr_tree: Tree containing all these subtrees

It creates one new internal node as a child of all_taxa_mrca_node
and places above mentioned subtrees as its children
"""
def MergeSubtrees(Curr_tree, first_cluster_mrca_node, second_cluster_mrca_node, all_taxa_mrca_node, taxa_list, Output_Text_File):
	if (DEBUG_LEVEL >= 2):
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
	Curr_tree.update_splits(delete_outdegree_one=False)

	if (DEBUG_LEVEL >= 2):
		fp = open(Output_Text_File, 'a')
		fp.write('\n label of newnode: ' + str(Node_Label(newnode)))
		fp.write('\n label of all taxa mrca node (recomputed): ' + str(Node_Label(Curr_tree.mrca(taxon_labels=taxa_list))))  
		fp.close()

	return Curr_tree

#-------------------------------------------
"""
this function merges two leaf nodes in the non-refined species tree
"""
def Merge_Leaves(Curr_tree, clust_species_list, idx1, idx2, taxa_list, Output_Text_File):
	"""
	both clusters are leaves - check whther there exist any level difference of them
	otherwise they can be treated as a sibling taxa pair in the species tree
	"""
	taxa1 = clust_species_list[idx1][0]
	taxa2 = clust_species_list[idx2][0]
	res = CheckR1RelationLevelBased(taxa1, taxa2)
	if (res == 0):
		Sibling_Couplet_List.append([taxa1, taxa2])

	clust1_mrca_node = Curr_tree.mrca(taxon_labels=clust_species_list[idx1])
	clust2_mrca_node = Curr_tree.mrca(taxon_labels=clust_species_list[idx2])
	all_taxa_mrca_node = Curr_tree.mrca(taxon_labels=taxa_list)
	Curr_tree = MergeSubtrees(Curr_tree, clust1_mrca_node, clust2_mrca_node, all_taxa_mrca_node, taxa_list, Output_Text_File)
	
	return Curr_tree

#-------------------------------------------
"""
this function merges one leaf node and another non leaf node (taxa cluster)
to refine the output species tree
"""
def Merge_Leaf_NonLeaf(Curr_tree, clust_species_list, leaf_idx, non_leaf_idx, taxa_list, Output_Text_File):
	
	# this is the MRCA node corresponding to the input leaf node
	leaf_taxon = clust_species_list[leaf_idx][0]
	leaf_mrca_node = Curr_tree.mrca(taxon_labels=clust_species_list[leaf_idx])
	
	# this is the MRCA node corresponding to the non leaf taxa cluster
	non_leaf_mrca_node = Curr_tree.mrca(taxon_labels=clust_species_list[non_leaf_idx])
	
	# MRCA node of the complete set of taxa (it can be original root of multifurcation as well)
	all_taxa_mrca_node = Curr_tree.mrca(taxon_labels=taxa_list)
	
	clust2_children = non_leaf_mrca_node.child_nodes()
	
	clust2_child1_taxa_list, clust2_child1_High_Level_Taxa_List, clust2_child1_Low_Level_Taxa_List = \
		FindLevelTaxa(Curr_tree, clust2_children[0], taxa_list)
	clust2_child2_taxa_list, clust2_child2_High_Level_Taxa_List, clust2_child2_Low_Level_Taxa_List = \
		FindLevelTaxa(Curr_tree, clust2_children[1], taxa_list)
	
	if (DEBUG_LEVEL >= 2):
		fp = open(Output_Text_File, 'a')
		fp.write('\n ===>>> Merging one leaf and other non leaf node') 
		fp.write('\n ---- First leaf idx list: ' + str(clust_species_list[leaf_idx]))
		fp.write('\n ---- Second non leaf idx list: ' + str(clust_species_list[non_leaf_idx]))
		fp.write('\n ---- clust2_child1_High_Level_Taxa_List: ' + str(clust2_child1_High_Level_Taxa_List)) 
		fp.write('\n ---- clust2_child1_Low_Level_Taxa_List: ' + str(clust2_child1_Low_Level_Taxa_List))
		fp.write('\n ---- clust2_child2_High_Level_Taxa_List: ' + str(clust2_child2_High_Level_Taxa_List)) 
		fp.write('\n ---- clust2_child2_Low_Level_Taxa_List: ' + str(clust2_child2_Low_Level_Taxa_List))
		fp.close()
	
	#----------------------------------------------------------------
	""" 
	representing first cluster as A (leaf node) and second cluster as (B,C)
	following operations can be performed only if B and C are not established sibling couplet
	"""
	#----------------------------------------------------------------
	if (CheckEstablishedSibling(clust2_child1_taxa_list, clust2_child2_taxa_list) == False):
		
		# trying the configuration (C, (A, B)) or (B, (A, C))
		res_clust2_child1_high_leaf_high = CheckR1RelationLevelBased(clust2_child1_High_Level_Taxa_List[0], leaf_taxon)
		val_clust2_child1_high_leaf_high = Get_R1R2Reln_ValDiff(clust2_child1_High_Level_Taxa_List[0], leaf_taxon)
		
		res_clust2_child2_high_leaf_high = CheckR1RelationLevelBased(clust2_child2_High_Level_Taxa_List[0], leaf_taxon)
		val_clust2_child2_high_leaf_high = Get_R1R2Reln_ValDiff(clust2_child2_High_Level_Taxa_List[0], leaf_taxon)
		
		if (DEBUG_LEVEL >= 2):
			fp = open(Output_Text_File, 'a')
			fp.write('\n ===>>> Obtaining pairwise results ---- ') 
			fp.write('\n ---- res_clust2_child1_high_leaf_high: ' + str(res_clust2_child1_high_leaf_high) + \
				' val_clust2_child1_high_leaf_high: ' + str(val_clust2_child1_high_leaf_high))
			fp.write('\n ---- res_clust2_child2_high_leaf_high: ' + str(res_clust2_child2_high_leaf_high) + \
				' val_clust2_child2_high_leaf_high: ' + str(val_clust2_child2_high_leaf_high))
			fp.close()
		
		# configuration (C, (A, B)) or (B, (A, C)) is present
		if (res_clust2_child1_high_leaf_high == 1) and (res_clust2_child2_high_leaf_high != 1):
			
			if (DEBUG_LEVEL >= 2):
				fp = open(Output_Text_File, 'a')
				fp.write('\n *** Condition A')
				fp.close()
			
			# (B, (A, C)) is present
			src_subtree_node = leaf_mrca_node
			dest_subtree_node = Curr_tree.mrca(taxon_labels=clust2_child2_taxa_list)
			Curr_tree = InsertSubTree(Curr_tree, src_subtree_node, dest_subtree_node, Output_Text_File)
			return Curr_tree
		
		if (res_clust2_child1_high_leaf_high != 1) and (res_clust2_child2_high_leaf_high == 1):
			
			if (DEBUG_LEVEL >= 2):
				fp = open(Output_Text_File, 'a')
				fp.write('\n *** Condition B')
				fp.close()
			
			# (C, (A, B)) is present
			src_subtree_node = leaf_mrca_node
			dest_subtree_node = Curr_tree.mrca(taxon_labels=clust2_child1_taxa_list)
			Curr_tree = InsertSubTree(Curr_tree, src_subtree_node, dest_subtree_node, Output_Text_File)
			return Curr_tree
		
		if (res_clust2_child1_high_leaf_high == 1) and (res_clust2_child2_high_leaf_high == 1):
			if (val_clust2_child1_high_leaf_high >= val_clust2_child2_high_leaf_high):
				
				if (DEBUG_LEVEL >= 2):
					fp = open(Output_Text_File, 'a')
					fp.write('\n *** Condition C')
					fp.close()
				
				# (B, (A, C)) is present
				src_subtree_node = leaf_mrca_node
				dest_subtree_node = Curr_tree.mrca(taxon_labels=clust2_child2_taxa_list)
				Curr_tree = InsertSubTree(Curr_tree, src_subtree_node, dest_subtree_node, Output_Text_File)
				return Curr_tree
			else:
				
				if (DEBUG_LEVEL >= 2):
					fp = open(Output_Text_File, 'a')
					fp.write('\n *** Condition D')
					fp.close()
				
				# (C, (A, B)) is present
				src_subtree_node = leaf_mrca_node
				dest_subtree_node = Curr_tree.mrca(taxon_labels=clust2_child1_taxa_list)
				Curr_tree = InsertSubTree(Curr_tree, src_subtree_node, dest_subtree_node, Output_Text_File)
				return Curr_tree
					
	# otherwise the configuration (A,(B,C)) will be employed
	Curr_tree = MergeSubtrees(Curr_tree, leaf_mrca_node, non_leaf_mrca_node, all_taxa_mrca_node, taxa_list, Output_Text_File)
	return Curr_tree
	
#-------------------------------------------
"""
this function has an input node and it traverses its underlying taxon and 
lists their level information
based on that, underlying taxon set is classified
"""
def FindLevelTaxa(Curr_tree, src_node, complete_taxa_list):
	taxa_list = []
	Level_List = []
	for n in src_node.leaf_nodes():
		if n.taxon.label in complete_taxa_list:
			taxa_list.append(n.taxon.label)
			Level_List.append(n.level())
			
	max_level = max(Level_List)
	min_level = min(Level_List)
	High_Level_Taxa_List = [taxa_list[i] for i, x in enumerate(Level_List) if x == min_level]
	Low_Level_Taxa_List = [taxa_list[i] for i, x in enumerate(Level_List) if x == max_level]
	
	#print 'initial High_Level_Taxa_List: ' + str(High_Level_Taxa_List)
	#print 'initial Low_Level_Taxa_List: ' + str(Low_Level_Taxa_List)
	""" 
	refine High_Level_Taxa_List and Low_Level_Taxa_List
	by applying tournament based selection to them
	"""
	while (len(High_Level_Taxa_List) > 1):
		del_list = []
		for i in range(0, len(High_Level_Taxa_List) - 1, 2):
			tax1 = High_Level_Taxa_List[i]
			tax2 = High_Level_Taxa_List[i+1]
			"""
			we check whether the couplet is already established as a sibling
			in such a case, we retain the taxa pair, without deleting it
			"""
			coup1 = [tax1, tax2]
			coup2 = [tax2, tax1]
			#print 'coup1: ' + str(coup1) + ' coup2: ' + str(coup2)
			if (coup1 not in Sibling_Couplet_List) and (coup2 not in Sibling_Couplet_List):
				res = CheckR1RelationLevelBased(tax1, tax2)
				if (res == 1):	# or (res == 0):
					del_list.append(tax2)
				elif (res == -1):
					del_list.append(tax1)
				# if res = 0 then both these taxa will be retained

		if (len(del_list) == 0):
			# nothing to delete - we can break from this loop
			# High_Level_Taxa_List[1:] = []
			break
		else:
			for d in del_list:
				High_Level_Taxa_List.remove(d)
		
	while (len(Low_Level_Taxa_List) > 1):
		del_list = []
		for i in range(0, len(Low_Level_Taxa_List) - 1, 2):
			tax1 = Low_Level_Taxa_List[i]
			tax2 = Low_Level_Taxa_List[i+1]
			"""
			we check whether the couplet is already established as a sibling
			in such a case, we retain the taxa pair, without deleting it
			"""
			coup1 = [tax1, tax2]
			coup2 = [tax2, tax1]
			#print 'coup1: ' + str(coup1) + ' coup2: ' + str(coup2)
			if (coup1 not in Sibling_Couplet_List) and (coup2 not in Sibling_Couplet_List):
				res = CheckR1RelationLevelBased(tax1, tax2)
				if (res == 1):	# or (res == 0):
					del_list.append(tax1)
				elif (res == -1):
					del_list.append(tax2)
				# if res = 0 then both these taxa will be retained

		if (len(del_list) == 0):
			# nothing to delete - we can break from this loop
			# Low_Level_Taxa_List[1:] = []
			break
		else:
			for d in del_list:
				Low_Level_Taxa_List.remove(d)
	
	#print 'final High_Level_Taxa_List: ' + str(High_Level_Taxa_List)
	#print 'final Low_Level_Taxa_List: ' + str(Low_Level_Taxa_List)
	
	return taxa_list, High_Level_Taxa_List, Low_Level_Taxa_List
	
#-------------------------------------------
"""
this function merges both non leaf nodes (taxa cluster) to refine the output species tree
"""
def Merge_Both_NonLeaf(Curr_tree, clust_species_list, idx1, idx2, taxa_list, Output_Text_File):
	
	# this is the MRCA node corresponding to the first cluster taxa list
	clust1_mrca_node = Curr_tree.mrca(taxon_labels=clust_species_list[idx1])
	
	# this is the MRCA node corresponding to the second cluster taxa list
	clust2_mrca_node = Curr_tree.mrca(taxon_labels=clust_species_list[idx2])
	
	# MRCA node of the complete set of taxa (it can be original root of multifurcation as well)
	all_taxa_mrca_node = Curr_tree.mrca(taxon_labels=taxa_list)
	
	# child nodes of both these MRCA nodes of respective taxa clusters
	clust1_children = clust1_mrca_node.child_nodes()
	clust1_child1_taxa_list, clust1_child1_High_Level_Taxa_List, clust1_child1_Low_Level_Taxa_List = \
		FindLevelTaxa(Curr_tree, clust1_children[0], taxa_list)
	clust1_child2_taxa_list, clust1_child2_High_Level_Taxa_List, clust1_child2_Low_Level_Taxa_List = \
		FindLevelTaxa(Curr_tree, clust1_children[1], taxa_list)
	
	if (DEBUG_LEVEL >= 2):
		fp = open(Output_Text_File, 'a')
		fp.write('\n ===>>> Merging both non leaf nodes') 
		fp.write('\n ---- First non leaf idx list: ' + str(clust_species_list[idx1]))
		fp.write('\n ---- clust1_child1_High_Level_Taxa_List: ' + str(clust1_child1_High_Level_Taxa_List)) 
		fp.write('\n ---- clust1_child1_Low_Level_Taxa_List: ' + str(clust1_child1_Low_Level_Taxa_List))
		fp.write('\n ---- clust1_child2_High_Level_Taxa_List: ' + str(clust1_child2_High_Level_Taxa_List)) 
		fp.write('\n ---- clust1_child2_Low_Level_Taxa_List: ' + str(clust1_child2_Low_Level_Taxa_List))
		fp.close()
	
	clust2_children = clust2_mrca_node.child_nodes()
	clust2_child1_taxa_list, clust2_child1_High_Level_Taxa_List, clust2_child1_Low_Level_Taxa_List = \
		FindLevelTaxa(Curr_tree, clust2_children[0], taxa_list)
	clust2_child2_taxa_list, clust2_child2_High_Level_Taxa_List, clust2_child2_Low_Level_Taxa_List = \
		FindLevelTaxa(Curr_tree, clust2_children[1], taxa_list)
	
	if (DEBUG_LEVEL >= 2):
		fp = open(Output_Text_File, 'a')
		fp.write('\n ===>>> Merging both non leaf nodes') 
		fp.write('\n ---- Second non leaf idx list: ' + str(clust_species_list[idx2]))
		fp.write('\n ---- clust2_child1_High_Level_Taxa_List: ' + str(clust2_child1_High_Level_Taxa_List)) 
		fp.write('\n ---- clust2_child1_Low_Level_Taxa_List: ' + str(clust2_child1_Low_Level_Taxa_List))
		fp.write('\n ---- clust2_child2_High_Level_Taxa_List: ' + str(clust2_child2_High_Level_Taxa_List)) 
		fp.write('\n ---- clust2_child2_Low_Level_Taxa_List: ' + str(clust2_child2_Low_Level_Taxa_List))
		fp.close()
	
	#----------------------------------------------------------------
	""" 
	representing first cluster as (A,B) and second cluster as (C,D)
	following operations can be performed only if A and B are not established sibling couplet
	"""
	#----------------------------------------------------------------
	if (CheckEstablishedSibling(clust1_child1_taxa_list, clust1_child2_taxa_list) == False):

		# trying the configuration (B, (A, (C, D))) or (A, (B, (C, D)))
		res_clust1_child1_low_clust2_child1_high = Check_R1Reln_Majority(clust1_child1_Low_Level_Taxa_List[0], clust2_child1_High_Level_Taxa_List[0])
		val_clust1_child1_low_clust2_child1_high = Get_R1R2Reln_ValDiff(clust1_child1_Low_Level_Taxa_List[0], clust2_child1_High_Level_Taxa_List[0])
		
		res_clust1_child1_low_clust2_child2_high = Check_R1Reln_Majority(clust1_child1_Low_Level_Taxa_List[0], clust2_child2_High_Level_Taxa_List[0])
		val_clust1_child1_low_clust2_child2_high = Get_R1R2Reln_ValDiff(clust1_child1_Low_Level_Taxa_List[0], clust2_child2_High_Level_Taxa_List[0])
		
		res_clust1_child2_low_clust2_child1_high  = Check_R1Reln_Majority(clust1_child2_Low_Level_Taxa_List[0], clust2_child1_High_Level_Taxa_List[0])
		val_clust1_child2_low_clust2_child1_high = Get_R1R2Reln_ValDiff(clust1_child2_Low_Level_Taxa_List[0], clust2_child1_High_Level_Taxa_List[0])
		
		res_clust1_child2_low_clust2_child2_high = Check_R1Reln_Majority(clust1_child2_Low_Level_Taxa_List[0], clust2_child2_High_Level_Taxa_List[0])
		val_clust1_child2_low_clust2_child2_high = Get_R1R2Reln_ValDiff(clust1_child2_Low_Level_Taxa_List[0], clust2_child2_High_Level_Taxa_List[0])
		
		if (DEBUG_LEVEL >= 2):
			fp = open(Output_Text_File, 'a')
			fp.write('\n ===>>> Obtaining pairwise results ---- ') 
			fp.write('\n ---- res_clust1_child1_low_clust2_child1_high: ' + str(res_clust1_child1_low_clust2_child1_high) + \
				' val_clust1_child1_low_clust2_child1_high: ' + str(val_clust1_child1_low_clust2_child1_high))
			fp.write('\n ---- res_clust1_child1_low_clust2_child2_high: ' + str(res_clust1_child1_low_clust2_child2_high) + \
				' val_clust1_child1_low_clust2_child2_high: ' + str(val_clust1_child1_low_clust2_child2_high))
			fp.write('\n ---- res_clust1_child2_low_clust2_child1_high: ' + str(res_clust1_child2_low_clust2_child1_high) + \
				' val_clust1_child2_low_clust2_child1_high: ' + str(val_clust1_child2_low_clust2_child1_high))
			fp.write('\n ---- res_clust1_child2_low_clust2_child2_high: ' + str(res_clust1_child2_low_clust2_child2_high) + \
				' val_clust1_child2_low_clust2_child2_high: ' + str(val_clust1_child2_low_clust2_child2_high))
			fp.close()
		
		# conditions in detail
		
		if (res_clust1_child1_low_clust2_child1_high == 1) and (res_clust1_child1_low_clust2_child2_high == 1) and \
			(res_clust1_child2_low_clust2_child1_high == 1) and (res_clust1_child2_low_clust2_child2_high == 1):
			
			# configuration (B, (A, (C, D))) or (A, (B, (C, D))) can be satisfied
				
			if (val_clust1_child1_low_clust2_child1_high < val_clust1_child2_low_clust2_child1_high) and \
				(val_clust1_child1_low_clust2_child2_high < val_clust1_child2_low_clust2_child2_high):
					
				if (DEBUG_LEVEL >= 2):
					fp = open(Output_Text_File, 'a')
					fp.write('\n *** Condition 1A')
					fp.close()
					
				src_subtree_node = clust2_mrca_node
				dest_subtree_node = Curr_tree.mrca(taxon_labels=clust1_child1_Low_Level_Taxa_List)
					
				# (B, (A, (C, D))) configuration is satisfied
				Curr_tree = InsertSubTree(Curr_tree, src_subtree_node, dest_subtree_node, Output_Text_File)
				
				return Curr_tree
			
			if (val_clust1_child1_low_clust2_child1_high > val_clust1_child2_low_clust2_child1_high) and \
				(val_clust1_child1_low_clust2_child2_high > val_clust1_child2_low_clust2_child2_high):

				if (DEBUG_LEVEL >= 2):
					fp = open(Output_Text_File, 'a')
					fp.write('\n *** Condition 1B')
					fp.close()

				src_subtree_node = clust2_mrca_node
				dest_subtree_node = Curr_tree.mrca(taxon_labels=clust1_child2_Low_Level_Taxa_List)

				# (A, (B, (C, D))) configuration is satisfied
				Curr_tree = InsertSubTree(Curr_tree, src_subtree_node, dest_subtree_node, Output_Text_File)
				
				return Curr_tree
	
	#----------------------------------------------------------------
	""" 
	representing first cluster as (A,B) and second cluster as (C,D)
	following operations can be performed only if C and D are not established sibling couplet
	"""
	#----------------------------------------------------------------
	if (CheckEstablishedSibling(clust2_child1_taxa_list, clust2_child2_taxa_list) == False):
	
		# trying the configuration (C, (D, (A, B))) or (D, (C, (A, B)))
		res_clust2_child1_low_clust1_child1_high = Check_R1Reln_Majority(clust2_child1_Low_Level_Taxa_List[0], clust1_child1_High_Level_Taxa_List[0])
		val_clust2_child1_low_clust1_child1_high = Get_R1R2Reln_ValDiff(clust2_child1_Low_Level_Taxa_List[0], clust1_child1_High_Level_Taxa_List[0])
		
		res_clust2_child1_low_clust1_child2_high = Check_R1Reln_Majority(clust2_child1_Low_Level_Taxa_List[0], clust1_child2_High_Level_Taxa_List[0])
		val_clust2_child1_low_clust1_child2_high = Get_R1R2Reln_ValDiff(clust2_child1_Low_Level_Taxa_List[0], clust1_child2_High_Level_Taxa_List[0])
		
		res_clust2_child2_low_clust1_child1_high = Check_R1Reln_Majority(clust2_child2_Low_Level_Taxa_List[0], clust1_child1_High_Level_Taxa_List[0])
		val_clust2_child2_low_clust1_child1_high = Get_R1R2Reln_ValDiff(clust2_child2_Low_Level_Taxa_List[0], clust1_child1_High_Level_Taxa_List[0])
		
		res_clust2_child2_low_clust1_child2_high = Check_R1Reln_Majority(clust2_child2_Low_Level_Taxa_List[0], clust1_child2_High_Level_Taxa_List[0])
		val_clust2_child2_low_clust1_child2_high = Get_R1R2Reln_ValDiff(clust2_child2_Low_Level_Taxa_List[0], clust1_child2_High_Level_Taxa_List[0])
		
		if (DEBUG_LEVEL >= 2):
			fp = open(Output_Text_File, 'a')
			fp.write('\n ===>>> Obtaining pairwise results ---- ') 
			fp.write('\n ---- res_clust2_child1_low_clust1_child1_high: ' + str(res_clust2_child1_low_clust1_child1_high) + \
				' val_clust2_child1_low_clust1_child1_high: ' + str(val_clust2_child1_low_clust1_child1_high))
			fp.write('\n ---- res_clust2_child1_low_clust1_child2_high: ' + str(res_clust2_child1_low_clust1_child2_high) + \
				' val_clust2_child1_low_clust1_child2_high: ' + str(val_clust2_child1_low_clust1_child2_high))
			fp.write('\n ---- res_clust2_child2_low_clust1_child1_high: ' + str(res_clust2_child2_low_clust1_child1_high) + \
				' val_clust2_child2_low_clust1_child1_high: ' + str(val_clust2_child2_low_clust1_child1_high))
			fp.write('\n ---- res_clust2_child2_low_clust1_child2_high: ' + str(res_clust2_child2_low_clust1_child2_high) + \
				' val_clust2_child2_low_clust1_child2_high: ' + str(val_clust2_child2_low_clust1_child2_high))
			fp.close()
		
		# conditions in detail 
		
		if (res_clust2_child1_low_clust1_child1_high == 1) and (res_clust2_child1_low_clust1_child2_high == 1) and \
			(res_clust2_child2_low_clust1_child1_high == 1) and (res_clust2_child2_low_clust1_child2_high == 1):
			
			# configuration (C, (D, (A, B))) or (D, (C, (A, B))) can be satisfied
				
			if (val_clust2_child1_low_clust1_child1_high < val_clust2_child2_low_clust1_child1_high) and \
				(val_clust2_child1_low_clust1_child2_high < val_clust2_child2_low_clust1_child2_high):
				
				if (DEBUG_LEVEL >= 2):
					fp = open(Output_Text_File, 'a')
					fp.write('\n *** Condition 2A')
					fp.close()
				
				src_subtree_node = clust1_mrca_node
				dest_subtree_node = Curr_tree.mrca(taxon_labels=clust2_child1_Low_Level_Taxa_List)
				
				# (D, (C, (A, B))) configuration is satisfied
				Curr_tree = InsertSubTree(Curr_tree, src_subtree_node, dest_subtree_node, Output_Text_File)
				
				return Curr_tree
			
			if (val_clust2_child1_low_clust1_child1_high > val_clust2_child2_low_clust1_child1_high) and \
				(val_clust2_child1_low_clust1_child2_high > val_clust2_child2_low_clust1_child2_high):

				if (DEBUG_LEVEL >= 2):
					fp = open(Output_Text_File, 'a')
					fp.write('\n *** Condition 2B')
					fp.close()

				src_subtree_node = clust1_mrca_node
				dest_subtree_node = Curr_tree.mrca(taxon_labels=clust2_child2_Low_Level_Taxa_List)

				# (C, (D, (A, B))) configuration is satisfied
				Curr_tree = InsertSubTree(Curr_tree, src_subtree_node, dest_subtree_node, Output_Text_File)
				
				return Curr_tree
	
	#-----------------------------
	#***************** add - sourya
	#-----------------------------
	#----------------------------------------------------------------
	""" 
	representing first cluster as (A,B) and second cluster as (C,D)
	following operations can be performed only if A and B are not established sibling couplet
	"""
	#----------------------------------------------------------------
	if (CheckEstablishedSibling(clust1_child1_taxa_list, clust1_child2_taxa_list) == False):
		
		# trying the configuration (B, (A, (C, D))) or (A, (B, (C, D)))
		
		res_clust1_child1_high_clust2_child1_high = CheckR1RelationLevelBased(clust1_child1_High_Level_Taxa_List[0], clust2_child1_High_Level_Taxa_List[0])
		val_clust1_child1_high_clust2_child1_high = Get_R1R2Reln_ValDiff(clust1_child1_High_Level_Taxa_List[0], clust2_child1_High_Level_Taxa_List[0])
		
		res_clust1_child1_high_clust2_child2_high = CheckR1RelationLevelBased(clust1_child1_High_Level_Taxa_List[0], clust2_child2_High_Level_Taxa_List[0])
		val_clust1_child1_high_clust2_child2_high = Get_R1R2Reln_ValDiff(clust1_child1_High_Level_Taxa_List[0], clust2_child2_High_Level_Taxa_List[0])
		
		res_clust1_child2_high_clust2_child1_high = CheckR1RelationLevelBased(clust1_child2_High_Level_Taxa_List[0], clust2_child1_High_Level_Taxa_List[0])
		val_clust1_child2_high_clust2_child1_high = Get_R1R2Reln_ValDiff(clust1_child2_High_Level_Taxa_List[0], clust2_child1_High_Level_Taxa_List[0])
		
		res_clust1_child2_high_clust2_child2_high = CheckR1RelationLevelBased(clust1_child2_High_Level_Taxa_List[0], clust2_child2_High_Level_Taxa_List[0])
		val_clust1_child2_high_clust2_child2_high = Get_R1R2Reln_ValDiff(clust1_child2_High_Level_Taxa_List[0], clust2_child2_High_Level_Taxa_List[0])
		
		if (DEBUG_LEVEL >= 2):
			fp = open(Output_Text_File, 'a')
			fp.write('\n ===>>> Obtaining pairwise results ---- ') 
			fp.write('\n ---- res_clust1_child1_high_clust2_child1_high: ' + str(res_clust1_child1_high_clust2_child1_high) + \
				' val_clust1_child1_high_clust2_child1_high: ' + str(val_clust1_child1_high_clust2_child1_high))
			fp.write('\n ---- res_clust1_child1_high_clust2_child2_high: ' + str(res_clust1_child1_high_clust2_child2_high) + \
				' val_clust1_child1_high_clust2_child2_high: ' + str(val_clust1_child1_high_clust2_child2_high))
			fp.write('\n ---- res_clust1_child2_high_clust2_child1_high: ' + str(res_clust1_child2_high_clust2_child1_high) + \
				' val_clust1_child2_high_clust2_child1_high: ' + str(val_clust1_child2_high_clust2_child1_high))
			fp.write('\n ---- res_clust1_child2_high_clust2_child2_high: ' + str(res_clust1_child2_high_clust2_child2_high) + \
				' val_clust1_child2_high_clust2_child2_high: ' + str(val_clust1_child2_high_clust2_child2_high))
			fp.close()
		
		# conditions in detail
		
		if (res_clust1_child1_high_clust2_child1_high == 1) and (res_clust1_child1_high_clust2_child2_high == 1) and \
			(res_clust1_child2_high_clust2_child1_high == 1) and (res_clust1_child2_high_clust2_child2_high == 1):
			
			# configuration (B, (A, (C, D))) or (A, (B, (C, D))) can be satisfied
				
			if (val_clust1_child1_high_clust2_child1_high >= val_clust1_child2_high_clust2_child1_high) and \
				(val_clust1_child1_high_clust2_child2_high >= val_clust1_child2_high_clust2_child2_high):
				
				"""
				here we satisfy the configuration (B, (A, (C, D)))
				note that > operator is used in condition checking
				also we use clust1_child1_taxa_list for MRCA computation
				"""
				
				if (DEBUG_LEVEL >= 2):
					fp = open(Output_Text_File, 'a')
					fp.write('\n *** Condition 3A')
					fp.close()
					
				src_subtree_node = clust2_mrca_node
				dest_subtree_node = Curr_tree.mrca(taxon_labels=clust1_child1_taxa_list)
					
				# (B, (A, (C, D))) configuration is satisfied
				Curr_tree = InsertSubTree(Curr_tree, src_subtree_node, dest_subtree_node, Output_Text_File)
				
				return Curr_tree
			
			if (val_clust1_child1_high_clust2_child1_high < val_clust1_child2_high_clust2_child1_high) and \
				(val_clust1_child1_high_clust2_child2_high < val_clust1_child2_high_clust2_child2_high):

				"""
				here we satisfy the configuration (A, (B, (C, D)))
				note that > operator is used in condition checking
				also we use clust1_child2_taxa_list for MRCA computation
				"""

				if (DEBUG_LEVEL >= 2):
					fp = open(Output_Text_File, 'a')
					fp.write('\n *** Condition 3B')
					fp.close()

				src_subtree_node = clust2_mrca_node
				dest_subtree_node = Curr_tree.mrca(taxon_labels=clust1_child2_taxa_list)

				# (A, (B, (C, D))) configuration is satisfied
				Curr_tree = InsertSubTree(Curr_tree, src_subtree_node, dest_subtree_node, Output_Text_File)
				
				return Curr_tree
	
	#----------------------------------------------------------------
	""" 
	representing first cluster as (A,B) and second cluster as (C,D)
	following operations can be performed only if C and D are not established sibling couplet
	"""
	#----------------------------------------------------------------
	if (CheckEstablishedSibling(clust2_child1_taxa_list, clust2_child2_taxa_list) == False):
		
		# trying the configuration (C, (D, (A, B))) or (D, (C, (A, B)))
		
		res_clust2_child1_high_clust1_child1_high = CheckR1RelationLevelBased(clust2_child1_High_Level_Taxa_List[0], clust1_child1_High_Level_Taxa_List[0])
		val_clust2_child1_high_clust1_child1_high = Get_R1R2Reln_ValDiff(clust2_child1_High_Level_Taxa_List[0], clust1_child1_High_Level_Taxa_List[0])
		
		res_clust2_child1_high_clust1_child2_high = CheckR1RelationLevelBased(clust2_child1_High_Level_Taxa_List[0], clust1_child2_High_Level_Taxa_List[0])
		val_clust2_child1_high_clust1_child2_high = Get_R1R2Reln_ValDiff(clust2_child1_High_Level_Taxa_List[0], clust1_child2_High_Level_Taxa_List[0])
		
		res_clust2_child2_high_clust1_child1_high = CheckR1RelationLevelBased(clust2_child2_High_Level_Taxa_List[0], clust1_child1_High_Level_Taxa_List[0])
		val_clust2_child2_high_clust1_child1_high = Get_R1R2Reln_ValDiff(clust2_child2_High_Level_Taxa_List[0], clust1_child1_High_Level_Taxa_List[0])
		
		res_clust2_child2_high_clust1_child2_high = CheckR1RelationLevelBased(clust2_child2_High_Level_Taxa_List[0], clust1_child2_High_Level_Taxa_List[0])
		val_clust2_child2_high_clust1_child2_high = Get_R1R2Reln_ValDiff(clust2_child2_High_Level_Taxa_List[0], clust1_child2_High_Level_Taxa_List[0])
		
		if (DEBUG_LEVEL >= 2):
			fp = open(Output_Text_File, 'a')
			fp.write('\n ===>>> Obtaining pairwise results ---- ') 
			fp.write('\n ---- res_clust2_child1_high_clust1_child1_high: ' + str(res_clust2_child1_high_clust1_child1_high) + \
				' val_clust2_child1_high_clust1_child1_high: ' + str(val_clust2_child1_high_clust1_child1_high))
			fp.write('\n ---- res_clust2_child1_high_clust1_child2_high: ' + str(res_clust2_child1_high_clust1_child2_high) + \
				' val_clust2_child1_high_clust1_child2_high: ' + str(val_clust2_child1_high_clust1_child2_high))
			fp.write('\n ---- res_clust2_child2_high_clust1_child1_high: ' + str(res_clust2_child2_high_clust1_child1_high) + \
				' val_clust2_child2_high_clust1_child1_high: ' + str(val_clust2_child2_high_clust1_child1_high))
			fp.write('\n ---- res_clust2_child2_high_clust1_child2_high: ' + str(res_clust2_child2_high_clust1_child2_high) + \
				' val_clust2_child2_high_clust1_child2_high: ' + str(val_clust2_child2_high_clust1_child2_high))
			fp.close()
		
		# conditions in detail 
		
		if (res_clust2_child1_high_clust1_child1_high == 1) and (res_clust2_child1_high_clust1_child2_high == 1) and \
			(res_clust2_child2_high_clust1_child1_high == 1) and (res_clust2_child2_high_clust1_child2_high == 1):
			
			# configuration (C, (D, (A, B))) or (D, (C, (A, B))) can be satisfied
				
			if (val_clust2_child1_high_clust1_child1_high >= val_clust2_child2_high_clust1_child1_high) and \
				(val_clust2_child1_high_clust1_child2_high >= val_clust2_child2_high_clust1_child2_high):
				
				"""
				here we satisfy the configuration (D, (C, (A, B)))
				note that > operator is used in condition checking
				also we use clust2_child1_taxa_list for MRCA computation
				"""
				
				if (DEBUG_LEVEL >= 2):
					fp = open(Output_Text_File, 'a')
					fp.write('\n *** Condition 4A')
					fp.close()
				
				src_subtree_node = clust1_mrca_node
				dest_subtree_node = Curr_tree.mrca(taxon_labels=clust2_child1_taxa_list)
				
				# (D, (C, (A, B))) configuration is satisfied
				Curr_tree = InsertSubTree(Curr_tree, src_subtree_node, dest_subtree_node, Output_Text_File)
				
				return Curr_tree
			
			if (val_clust2_child1_high_clust1_child1_high < val_clust2_child2_high_clust1_child1_high) and \
				(val_clust2_child1_high_clust1_child2_high < val_clust2_child2_high_clust1_child2_high):

				"""
				here we satisfy the configuration (C, (D, (A, B)))
				note that > operator is used in condition checking
				also we use clust2_child2_taxa_list for MRCA computation
				"""

				if (DEBUG_LEVEL >= 2):
					fp = open(Output_Text_File, 'a')
					fp.write('\n *** Condition 4B')
					fp.close()

				src_subtree_node = clust1_mrca_node
				dest_subtree_node = Curr_tree.mrca(taxon_labels=clust2_child2_taxa_list)

				# (C, (D, (A, B))) configuration is satisfied
				Curr_tree = InsertSubTree(Curr_tree, src_subtree_node, dest_subtree_node, Output_Text_File)
				
				return Curr_tree
	
	#-----------------------------
	#***************** end add - sourya
	#-----------------------------
	
	# otherwise the configuration ((A,B),(C,D)) will be employed
	Curr_tree = MergeSubtrees(Curr_tree, clust1_mrca_node, clust2_mrca_node, all_taxa_mrca_node, taxa_list, Output_Text_File)
	
	return Curr_tree	

#-------------------------------------------
"""
this is a new function to merge taxa clusters for agglomerative clustering - sourya
"""
def Merge_Cluster_Pair_New(Curr_tree, clust_species_list, min_idx_i, min_idx_j, taxa_list, Output_Text_File):

	isleaf_clust1 = IsLeafCluster(clust_species_list, min_idx_i)
	isleaf_clust2 = IsLeafCluster(clust_species_list, min_idx_j)
	if (isleaf_clust1 == True) and (isleaf_clust2 == True):
		Curr_tree = Merge_Leaves(Curr_tree, clust_species_list, min_idx_i, min_idx_j, taxa_list, Output_Text_File)
	elif (isleaf_clust1 == True) and (isleaf_clust2 == False):
		Curr_tree = Merge_Leaf_NonLeaf(Curr_tree, clust_species_list, min_idx_i, min_idx_j, taxa_list, Output_Text_File)
	elif (isleaf_clust1 == False) and (isleaf_clust2 == True):
		Curr_tree = Merge_Leaf_NonLeaf(Curr_tree, clust_species_list, min_idx_j, min_idx_i, taxa_list, Output_Text_File)
	elif (isleaf_clust1 == False) and (isleaf_clust2 == False):
		Curr_tree = Merge_Both_NonLeaf(Curr_tree, clust_species_list, min_idx_i, min_idx_j, taxa_list, Output_Text_File)
	
	return Curr_tree

#-------------------------------------------
"""
this function merges a pair of clusters whose indices are pointed by the min_idx_i and min_idx_j entries
this is part of the proposed agglomerative clustering
taxa_list is the union of these two clusters (species contents)
"""
def Merge_Cluster_Pair(Curr_tree, clust_species_list, min_idx_i, min_idx_j, taxa_list, Output_Text_File):
	isleaf_clust1 = IsLeafCluster(clust_species_list, min_idx_i)
	isleaf_clust2 = IsLeafCluster(clust_species_list, min_idx_j)
	
	if (isleaf_clust1):
		first_cluster_mrca_node = Curr_tree.find_node_with_taxon_label(clust_species_list[min_idx_i][0])
	else:
		first_cluster_mrca_node = Curr_tree.mrca(taxon_labels=clust_species_list[min_idx_i])
	
	if (isleaf_clust2):
		second_cluster_mrca_node = Curr_tree.find_node_with_taxon_label(clust_species_list[min_idx_j][0])
	else:
		second_cluster_mrca_node = Curr_tree.mrca(taxon_labels=clust_species_list[min_idx_j])
	
	all_taxa_mrca_node = Curr_tree.mrca(taxon_labels=taxa_list)
	
	Curr_tree = MergeSubtrees(Curr_tree, first_cluster_mrca_node, second_cluster_mrca_node, all_taxa_mrca_node, taxa_list, Output_Text_File)
	
	return Curr_tree

#-------------------------------------------
"""
this function processes one internal node (basically the children list)
to resolve multifurcation
"""
def ResolveMultifurcation(Curr_tree, clust_species_list, no_of_input_clusters, Output_Text_File, NJ_RULE_USED):
	# total number of clusters
	no_of_clust = no_of_input_clusters

	# allocate a 2D square matrix of no_of_clust dimension
	DistMat = numpy.zeros((no_of_clust, no_of_clust), dtype=numpy.float)
	Norm_DistMat = numpy.zeros((no_of_clust, no_of_clust), dtype=numpy.float)

	# here we compute the ILS score of the current cluster pair
	# with respect to input gene tree list
	if (DEBUG_LEVEL > 2):
		fp = open(Output_Text_File, 'a')
		fp.write('\n Examining ILS score for individual cluster pairs ')
		fp.close()      

	#---------------------------------------
	# using single taxon as a representative of the taxa cluster
	#Fill_DistMat_SingleEntry(DistMat, no_of_clust, clust_species_list)
	
	# using average information from a taxa cluster - currently commented - requires more running time
	Fill_DistMat_AvgEntry(DistMat, no_of_clust, clust_species_list)
	#---------------------------------------

	# loop to execute the agglomerative clustering
	while(no_of_clust > 2):
		min_idx_i, min_idx_j = Get_NJ_Based_Min_Pair_Idx(DistMat, Norm_DistMat, \
			no_of_clust, clust_species_list, NJ_RULE_USED, Output_Text_File)

		# note down the taxa list in these two indices of the clust_species_list
		taxa_list = []
		for x in clust_species_list[min_idx_i]:
			taxa_list.append(x)
		for x in clust_species_list[min_idx_j]:
			taxa_list.append(x)

		if (DEBUG_LEVEL >= 2):
			fp = open(Output_Text_File, 'a')
			fp.write('\n min_idx_i ' + str(min_idx_i) + ' min_idx_j : ' + str(min_idx_j))
			fp.write('\n min_idx_i species list ' + str(clust_species_list[min_idx_i]))
			fp.write('\n min_idx_j species list ' + str(clust_species_list[min_idx_j]))
			fp.write('\n complete species list ' + str(taxa_list))
			fp.close()
			
		"""
		now we merge the pair of clusters pointed by these indices
		"""
		Curr_tree = Merge_Cluster_Pair(Curr_tree, clust_species_list, min_idx_i, min_idx_j, taxa_list, Output_Text_File)
		#Curr_tree = Merge_Cluster_Pair_New(Curr_tree, clust_species_list, min_idx_i, min_idx_j, taxa_list, Output_Text_File)
		#---------------------------------------------------------      
		"""
		adjust the DistMat by inserting one new row and column corresponding to the new cluster
		and then deleting the information of earlier two clusters
		"""
		# first append one row
		DistMat = numpy.vstack((DistMat, numpy.zeros((1, no_of_clust), dtype=numpy.float)))
		# then append one column
		DistMat = numpy.hstack((DistMat, numpy.zeros((no_of_clust + 1, 1), dtype=numpy.float)))
		# now apply reshape operation to get proper square matrix dimension
		DistMat = numpy.reshape(DistMat, ((no_of_clust + 1), (no_of_clust + 1)), order='C')
		
		# now fill the elements of the new added row and column
		for k in range(no_of_clust):
			# for any index k, the number of extra lineage is the maximum of these three quantities
			if (NJ_RULE_USED == AGGLO_CLUST):
				DistMat[k][no_of_clust] = max(DistMat[k][min_idx_i], DistMat[k][min_idx_j], DistMat[min_idx_i][min_idx_j])
				#DistMat[k][no_of_clust] = (DistMat[k][min_idx_i] * len(clust_species_list[min_idx_i]) + DistMat[k][min_idx_j] * len(clust_species_list[min_idx_j])) / (len(clust_species_list[min_idx_i]) + len(clust_species_list[min_idx_j]))
				#DistMat[k][no_of_clust] = max(DistMat[k][min_idx_i], DistMat[k][min_idx_j])	# modified - sourya
			else:
				DistMat[k][no_of_clust] = (DistMat[k][min_idx_i] + DistMat[k][min_idx_j] - DistMat[min_idx_i][min_idx_j]) / 2.0
			# symmetric property
			DistMat[no_of_clust][k] = DistMat[k][no_of_clust]
		
		# now remove the rows and columns corresponding to min_idx_i and min_idx_j
		DistMat = numpy.delete(DistMat, (min_idx_i), axis=0)	# delete the row
		DistMat = numpy.delete(DistMat, (min_idx_i), axis=1)	# delete the column
		DistMat = numpy.delete(DistMat, (min_idx_j - 1), axis=0)	# delete the row
		DistMat = numpy.delete(DistMat, (min_idx_j - 1), axis=1)	# delete the column

		# clear Norm_DistMat
		Norm_DistMat = numpy.delete(Norm_DistMat, (min_idx_i), axis=0)	# delete the row
		Norm_DistMat = numpy.delete(Norm_DistMat, (min_idx_i), axis=1)	# delete the column
		Norm_DistMat.fill(0)
		
		# remove individual clusters' taxa information from the clust_species_list
		# and add taxa_list as a new element
		clust_species_list.pop(min_idx_i)
		clust_species_list.pop(min_idx_j - 1)
		clust_species_list.append(taxa_list)    
		
		# decrement the number of clusters considered
		no_of_clust = no_of_clust - 1
	
	return
            
#-------------------------------------------
# this function refines input supertree such that the supertree becomes binary
# this is required for proper benchmarking with existing binary tree construction methods on 
# ILS sorting
def Refine_Supertree_Binary_Form(Curr_tree, Output_Text_File, NJ_RULE_USED):

	# add - sourya
	# maintain the list of couplets
	for curr_node in Curr_tree.postorder_internal_node_iter():
		curr_node_children = curr_node.child_nodes()
		if (len(curr_node_children) == 2) and (curr_node_children[0].is_leaf() == True) and (curr_node_children[1].is_leaf() == True):
			subl = []
			for x in curr_node_children:
				if (x.is_leaf() == True):
					subl.append(x.taxon.label)
			Sibling_Couplet_List.append(subl)
	# end add - sourya
	
	if (DEBUG_LEVEL >= 2):
		fp = open(Output_Text_File, 'a')
		fp.write('\n Sibling_Couplet_List ' + str(Sibling_Couplet_List))
		fp.close()
	
	"""
	we traverse input tree internal nodes in postorder fashion
	and list the child nodes of it
	if the no of children > 2 then it is a case of multifurcation
	for resolving
	"""
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
			ResolveMultifurcation(Curr_tree, clust_species_list, len(curr_node_children), Output_Text_File, NJ_RULE_USED)
