#!/usr/bin/env python

import Header
from Header import *
import UtilFunc
from UtilFunc import *

#-------------------------------------------
"""
this function fills the distance matrix using normalized excess gene count
for a particular taxa cluster, it uses one representative taxon of that taxa cluster
@param type_of_output: if 0, computes the average of XL measures
													1, returns the minimum of XL measures
													2, returns the maximum of XL measures
"""
def Fill_DistMat_SingleEntry(DistMat, no_of_clust, clust_species_list, DIST_MAT_TYPE, type_of_output=0):
	"""
	check and explore each pair of taxa clusters
	"""
	for i in range(no_of_clust - 1):
		#"""
		#arrange elements of clust_species_list[i] 
		#in preorder manner
		#"""
		## this is the MRCA node corresponding to the clust_species_list[i]
		#clust_i_mrca_node = Curr_tree.mrca(taxon_labels=clust_species_list[i])
		#clust_i_preorder_taxa_list = []
		#for n in clust_i_mrca_node.preorder_iter():
			#if (n.is_leaf() == True):
				#if n.taxon.label in clust_species_list[i]:
					#clust_i_preorder_taxa_list.append(n.taxon.label)
		
		for j in range(i+1, no_of_clust):
			#"""
			#clust_species_list[i] and clust_species_list[j]
			#contain two taxa list of one or more elements
			#"""
			#""" 
			#first we arrange elements of clust_species_list[i] and clust_species_list[j] 
			#in preorder manner
			#"""
			## this is the MRCA node corresponding to the clust_species_list[j]
			#clust_j_mrca_node = Curr_tree.mrca(taxon_labels=clust_species_list[j])
	
			
			#clust_j_preorder_taxa_list = []
	
	
			#for n in clust_j_mrca_node.preorder_iter():
				#if (n.is_leaf() == True):
					#if n.taxon.label in clust_species_list[j]:
						#clust_j_preorder_taxa_list.append(n.taxon.label)

			entry = FindAvgXL(clust_species_list[i], clust_species_list[j], DIST_MAT_TYPE, True, type_of_output)
			DistMat[j][i] = DistMat[i][j] = entry

	return
	
#-------------------------------------------
"""
this function fills the distance matrix using normalized excess gene count
for a particular taxa cluster, it uses aerage information of that taxa cluster
@param type_of_output: if 0, computes the average of XL measures
													1, returns the minimum of XL measures
													2, returns the maximum of XL measures
"""
def Fill_DistMat_AvgEntry(DistMat, no_of_clust, clust_species_list, DIST_MAT_TYPE, type_of_output=0):
	"""
	check and explore each pair of taxa clusters
	"""
	for i in range(no_of_clust - 1):
		for j in range(i+1, no_of_clust):
			# here both clust_species_list[i] and clust_species_list[j]
			# are one element lists (according to their construction)
			# we have extracted the corresponding element by using [0] operator (extracting first element)
			entry = FindAvgXL(clust_species_list[i], clust_species_list[j], DIST_MAT_TYPE, False, type_of_output)
			DistMat[j][i] = DistMat[i][j] = entry
	
	return

#-------------------------------------------
#"""
#following function checks if we are merging two leaves 
#which do not have R3 as their consensus relation
#in such a case, the function returns FALSE
#"""
#def CheckMergeImproperLeafandNonleaf(clust_species_list, i, j):
	#if (IsLeafCluster(clust_species_list, i) == True) and (IsLeafCluster(clust_species_list, j) == False):
		#leaf_idx = i
		#non_leaf_idx = j
	#elif (IsLeafCluster(clust_species_list, i) == False) and (IsLeafCluster(clust_species_list, j) == True):
		#leaf_idx = j
		#non_leaf_idx = i
	#else:
		#return False
	
	#leaf_taxon = clust_species_list[leaf_idx][0]
	#for nonleaf_taxon in clust_species_list[non_leaf_idx]:
		#key1 = (leaf_taxon, nonleaf_taxon)
		#key2 = (nonleaf_taxon, leaf_taxon)
		#if key1 in TaxaPair_Reln_Dict:
			#if (TaxaPair_Reln_Dict[key1]._GetR1R2LevelDiff(False) <= 0) or \
				#(TaxaPair_Reln_Dict[key1]._GetRelnLevelDiff(RELATION_R1, RELATION_R2, False, False) <= 0):
				#return True
		#if key2 in TaxaPair_Reln_Dict:
			#if (TaxaPair_Reln_Dict[key2]._GetR1R2LevelDiff(False) >= 0) or \
				#(TaxaPair_Reln_Dict[key2]._GetRelnLevelDiff(RELATION_R1, RELATION_R2, False, False) >= 0):
				#return True

	#return False
	
##-------------------------------------------
#"""
#following function checks if we are merging two leaves 
#which do not have R3 as their consensus relation
#in such a case, the function returns FALSE
#"""
#def CheckMergeImproperLeaves(clust_species_list, i, j):
	#if (IsLeafCluster(clust_species_list, i) == True) and (IsLeafCluster(clust_species_list, j) == True):
		#key1 = (clust_species_list[i][0], clust_species_list[j][0])
		#key2 = (clust_species_list[j][0], clust_species_list[i][0])
		#if key1 in TaxaPair_Reln_Dict:
			#if (TaxaPair_Reln_Dict[key1]._CheckTargetRelnConsensus(RELATION_R3) == False):
				#return True
		#if key2 in TaxaPair_Reln_Dict:
			#if (TaxaPair_Reln_Dict[key2]._CheckTargetRelnConsensus(RELATION_R3) == False):
				#return True
	
	#return False

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
				if 1:	#(CheckMergeImproperLeaves(clust_species_list, i, j) == False) and (CheckMergeImproperLeafandNonleaf(clust_species_list, i, j) == False):	# add - sourya
					min_val = Norm_DistMat[i][j]
					min_idx_i = i
					min_idx_j = j
			elif (Norm_DistMat[i][j] == min_val):
				if 1:	#(CheckMergeImproperLeaves(clust_species_list, i, j) == False) and (CheckMergeImproperLeafandNonleaf(clust_species_list, i, j) == False):	# add - sourya
					# comment - sourya
					## here we prioritize the cluster pair having minimum number of species
					#if (len(clust_species_list[i]) + len(clust_species_list[j])) < (len(clust_species_list[min_idx_i]) + len(clust_species_list[min_idx_j])):
						#min_idx_i = i
						#min_idx_j = j
					# end comment - sourya
					# add - sourya
					#if (Get_R1R2Reln_ValDiff(clust_species_list[i][0], clust_species_list[j][0], True) < \
						#Get_R1R2Reln_ValDiff(clust_species_list[min_idx_i][0], clust_species_list[min_idx_j][0], True)):
					if (GetR3RelnLevelConsVal(clust_species_list[i][0], clust_species_list[j][0]) > \
						GetR3RelnLevelConsVal(clust_species_list[min_idx_i][0], clust_species_list[min_idx_j][0])):
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


##-------------------------------------------------------------
#"""
#this function checks whether the input pair of taxa list is basically 
#a couplet, belonging to an existing set of siblings
#"""
#def CheckEstablishedSibling(taxa_list1, taxa_list2):
	#if (len(taxa_list1) == 1) and (len(taxa_list2) == 1):
		#coup1 = [taxa_list1[0], taxa_list2[0]]
		#coup2 = [taxa_list2[0], taxa_list1[0]]
		##print 'coup1: ' + str(coup1) + '  or coup2: ' + str(coup2) + ' is an established sibling pair'
		#if (coup1 in Sibling_Couplet_List) or (coup2 in Sibling_Couplet_List):
			#return True
	#return False
	
##-----------------------------------------------------------------
#def Get_R1R2Reln_ValDiff(x1, x2, abs_comp):
	#key1 = (x1, x2)
	#key2 = (x2, x1)
	#if key1 in TaxaPair_Reln_Dict:	
		#return TaxaPair_Reln_Dict[key1]._GetR1R2LevelDiff(abs_comp)
	#elif key2 in TaxaPair_Reln_Dict:
		#if (abs_comp == False):
			#return (-1) * TaxaPair_Reln_Dict[key2]._GetR1R2LevelDiff(abs_comp)
		#else:
			#return TaxaPair_Reln_Dict[key2]._GetR1R2LevelDiff(abs_comp)
	
	#return 0

#-----------------------------------------------------------------
def GetR3RelnLevelConsVal(x1, x2):
	key1 = (x1, x2)
	key2 = (x2, x1)
	if key1 in TaxaPair_Reln_Dict:	
		return TaxaPair_Reln_Dict[key1]._GetRelnLevelDiff(RELATION_R3, -1, False, True)
	elif key2 in TaxaPair_Reln_Dict:
		return TaxaPair_Reln_Dict[key2]._GetRelnLevelDiff(RELATION_R3, -1, False, True)

##-----------------------------------------------------------------
#"""
#this function checks a couplet whether they can be related by R1 relation
#according to the all relation level count difference statistic
#"""
#def Check_R1Reln_Majority(x1, x2):
	#key1 = (x1, x2)
	#key2 = (x2, x1)
	#if key1 in TaxaPair_Reln_Dict:
		#if (TaxaPair_Reln_Dict[key1]._CheckTargetRelnLevelConsensus(RELATION_R2)):
			#return -1
		#elif (TaxaPair_Reln_Dict[key1]._CheckTargetRelnLevelConsensus(RELATION_R1)):
			#return 1
		#elif (TaxaPair_Reln_Dict[key1]._CheckHigherTargetRelnLevelValue(RELATION_R2)):
			#return -2
		#elif (TaxaPair_Reln_Dict[key1]._CheckHigherTargetRelnLevelValue(RELATION_R1)):
			#return 2
	#elif key2 in TaxaPair_Reln_Dict:
		#if (TaxaPair_Reln_Dict[key2]._CheckTargetRelnLevelConsensus(RELATION_R2)):
			#return 1
		#elif (TaxaPair_Reln_Dict[key2]._CheckTargetRelnLevelConsensus(RELATION_R1)):
			#return -1
		#if (TaxaPair_Reln_Dict[key2]._CheckHigherTargetRelnLevelValue(RELATION_R2)):
			#return 2
		#elif (TaxaPair_Reln_Dict[key2]._CheckHigherTargetRelnLevelValue(RELATION_R1)):
			#return -2
	
	#return KEY_ABSENCE_INDICATOR	# key absence indicator
	
##-----------------------------------------------------------------
#"""
#this function checks according to the all relation level count difference statistic
#which cluster can be ancestor and which can be descendant
#returns - 1) binary value 1 or 0 depending on the couplet is a sibling or not
#2) indices x and y where the species in index x is placed at a higher level than the species at index y
#with respect to the input gene trees' configurations
#"""

## this function is modified by - sourya - 29.01.2016

#def CheckR1RelationLevelBased(x1, x2):
	#key1 = (x1, x2)
	#key2 = (x2, x1)
	#if key1 in TaxaPair_Reln_Dict:
		#if (TaxaPair_Reln_Dict[key1]._CheckTargetRelnLevelConsensus(RELATION_R3)):
			#return 0
		##elif (TaxaPair_Reln_Dict[key1]._CheckHigherTargetRelnLevelValue(RELATION_R2)):
			##return -1
		##elif (TaxaPair_Reln_Dict[key1]._CheckHigherTargetRelnLevelValue(RELATION_R1)):
			##return 1
		#elif (TaxaPair_Reln_Dict[key1]._CheckTargetRelnLevelConsensus(RELATION_R2)):
			#return -1
		#elif (TaxaPair_Reln_Dict[key1]._CheckTargetRelnLevelConsensus(RELATION_R1)):
			#return 1
	#elif key2 in TaxaPair_Reln_Dict:
		#if (TaxaPair_Reln_Dict[key2]._CheckTargetRelnLevelConsensus(RELATION_R3)):
			#return 0
		##elif (TaxaPair_Reln_Dict[key2]._CheckHigherTargetRelnLevelValue(RELATION_R2)):
			##return 1
		##elif (TaxaPair_Reln_Dict[key2]._CheckHigherTargetRelnLevelValue(RELATION_R1)):
			##return -1
		#elif (TaxaPair_Reln_Dict[key2]._CheckTargetRelnLevelConsensus(RELATION_R2)):
			#return 1
		#elif (TaxaPair_Reln_Dict[key2]._CheckTargetRelnLevelConsensus(RELATION_R1)):
			#return -1
	
	#return KEY_ABSENCE_INDICATOR	# key absence indicator

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
	#res = CheckR1RelationLevelBased(taxa1, taxa2)
	#if (res == 0):
		#Sibling_Couplet_List.append([taxa1, taxa2])

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
def Merge_Leaf_NonLeaf(Curr_tree, clust_species_list, leaf_idx, non_leaf_idx, taxa_list, Output_Text_File, DIST_MAT_TYPE, DIST_MAT_UPDATE):
	
	#----------------------------------------------------------------
	# representing first cluster as A (leaf node) and second cluster as (B,C)
	#----------------------------------------------------------------
	
	# this is the MRCA node corresponding to the input leaf node
	leaf_taxon = clust_species_list[leaf_idx][0]
	leaf_mrca_node = Curr_tree.mrca(taxon_labels=clust_species_list[leaf_idx])
	
	# this is the MRCA node corresponding to the non leaf taxa cluster
	non_leaf_mrca_node = Curr_tree.mrca(taxon_labels=clust_species_list[non_leaf_idx])
	
	# MRCA node of the complete set of taxa (it can be original root of multifurcation as well)
	all_taxa_mrca_node = Curr_tree.mrca(taxon_labels=taxa_list)
	
	clust2_children = non_leaf_mrca_node.child_nodes()
	clust2_child1_taxa_list = []
	clust2_child2_taxa_list = []
	
	# comment - sourya
	#for n in clust2_children[0].leaf_nodes():
		#if n.taxon.label in taxa_list:
			#clust2_child1_taxa_list.append(n.taxon.label)
	# add - sourya
	for n in clust2_children[0].preorder_iter():
		if (n.is_leaf() == True):
			if n.taxon.label in taxa_list:
				clust2_child1_taxa_list.append(n.taxon.label)
	# comment - sourya
	#for n in clust2_children[1].leaf_nodes():
		#if n.taxon.label in taxa_list:
			#clust2_child2_taxa_list.append(n.taxon.label)
	# add - sourya
	for n in clust2_children[1].preorder_iter():
		if (n.is_leaf() == True):
			if n.taxon.label in taxa_list:
				clust2_child2_taxa_list.append(n.taxon.label)
	
	if (DEBUG_LEVEL >= 2):
		fp = open(Output_Text_File, 'a')
		fp.write('\n ===>>> Merging one leaf and other non leaf node') 
		fp.write('\n ---- First leaf idx list: ' + str(clust_species_list[leaf_idx]))
		fp.write('\n ---- Second non leaf idx list: ' + str(clust_species_list[non_leaf_idx]))
		fp.write('\n ---- clust2_child1_taxa_list: ' + str(clust2_child1_taxa_list)) 
		fp.write('\n ---- clust2_child2_taxa_list: ' + str(clust2_child2_taxa_list))
		fp.close()
	
	# obtain three sets of XL values
	XL_clust2_child1_clust2_child2 = FindAvgXL(clust2_child1_taxa_list, clust2_child2_taxa_list, DIST_MAT_TYPE, False, (DIST_MAT_UPDATE - 1))
	XL_clust2_child1_leaf = FindAvgXL(clust2_child1_taxa_list, clust_species_list[leaf_idx], DIST_MAT_TYPE, False, (DIST_MAT_UPDATE - 1))
	XL_clust2_child2_leaf = FindAvgXL(clust2_child2_taxa_list, clust_species_list[leaf_idx], DIST_MAT_TYPE, False, (DIST_MAT_UPDATE - 1))

	if (DEBUG_LEVEL >= 2):
		fp = open(Output_Text_File, 'a')
		fp.write('\n ---- XL_clust2_child1_clust2_child2: ' + str(XL_clust2_child1_clust2_child2))
		fp.write('\n ---- XL_clust2_child1_leaf: ' + str(XL_clust2_child1_leaf))
		fp.write('\n ---- XL_clust2_child2_leaf: ' + str(XL_clust2_child2_leaf))
		fp.close()
		
	## add - sourya
	#"""
	#if one of the children of the clust2 is a leaf
	#then we proceed for the simple merging
	#"""
	#if (len(clust2_child1_taxa_list) == 1) or (len(clust2_child2_taxa_list) == 1):
		## configuration (A,(B,C)) will be employed
		#Curr_tree = MergeSubtrees(Curr_tree, leaf_mrca_node, non_leaf_mrca_node, all_taxa_mrca_node, taxa_list, Output_Text_File)
		#return Curr_tree
	## end add - sourya

	XL_list = [XL_clust2_child1_clust2_child2, XL_clust2_child1_leaf, XL_clust2_child2_leaf]
	if (XL_clust2_child1_clust2_child2 == min(XL_list)):
		# configuration (A,(B,C)) will be employed
		Curr_tree = MergeSubtrees(Curr_tree, leaf_mrca_node, non_leaf_mrca_node, all_taxa_mrca_node, taxa_list, Output_Text_File)
	elif (XL_clust2_child1_leaf == min(XL_list)):
		# configuration (B, (A, C)) will be employed
		src_subtree_node = leaf_mrca_node
		dest_subtree_node = Curr_tree.mrca(taxon_labels=clust2_child2_taxa_list)
		Curr_tree = InsertSubTree(Curr_tree, src_subtree_node, dest_subtree_node, Output_Text_File)
	else:
		# configuration (C, (A, B)) will be employed 
		src_subtree_node = leaf_mrca_node
		dest_subtree_node = Curr_tree.mrca(taxon_labels=clust2_child1_taxa_list)
		Curr_tree = InsertSubTree(Curr_tree, src_subtree_node, dest_subtree_node, Output_Text_File)
	
	return Curr_tree

#-------------------------------------------
#"""
#this function merges one leaf node and another non leaf node (taxa cluster)
#to refine the output species tree
#"""
#def Merge_Leaf_NonLeaf_OLD(Curr_tree, clust_species_list, leaf_idx, non_leaf_idx, taxa_list, Output_Text_File):
	
	## this is the MRCA node corresponding to the input leaf node
	#leaf_taxon = clust_species_list[leaf_idx][0]
	#leaf_mrca_node = Curr_tree.mrca(taxon_labels=clust_species_list[leaf_idx])
	
	## this is the MRCA node corresponding to the non leaf taxa cluster
	#non_leaf_mrca_node = Curr_tree.mrca(taxon_labels=clust_species_list[non_leaf_idx])
	
	## MRCA node of the complete set of taxa (it can be original root of multifurcation as well)
	#all_taxa_mrca_node = Curr_tree.mrca(taxon_labels=taxa_list)
	
	#clust2_children = non_leaf_mrca_node.child_nodes()
	
	#clust2_child1_taxa_list, clust2_child1_High_Level_Taxa_List, clust2_child1_Low_Level_Taxa_List = \
		#FindLevelTaxa(Curr_tree, clust2_children[0], taxa_list)
	#clust2_child2_taxa_list, clust2_child2_High_Level_Taxa_List, clust2_child2_Low_Level_Taxa_List = \
		#FindLevelTaxa(Curr_tree, clust2_children[1], taxa_list)
	
	#if (DEBUG_LEVEL >= 2):
		#fp = open(Output_Text_File, 'a')
		#fp.write('\n ===>>> Merging one leaf and other non leaf node') 
		#fp.write('\n ---- First leaf idx list: ' + str(clust_species_list[leaf_idx]))
		#fp.write('\n ---- Second non leaf idx list: ' + str(clust_species_list[non_leaf_idx]))
		#fp.write('\n ---- clust2_child1_High_Level_Taxa_List: ' + str(clust2_child1_High_Level_Taxa_List)) 
		#fp.write('\n ---- clust2_child1_Low_Level_Taxa_List: ' + str(clust2_child1_Low_Level_Taxa_List))
		#fp.write('\n ---- clust2_child2_High_Level_Taxa_List: ' + str(clust2_child2_High_Level_Taxa_List)) 
		#fp.write('\n ---- clust2_child2_Low_Level_Taxa_List: ' + str(clust2_child2_Low_Level_Taxa_List))
		#fp.close()
	
	##----------------------------------------------------------------
	#""" 
	#representing first cluster as A (leaf node) and second cluster as (B,C)
	#following operations can be performed only if B and C are not established sibling couplet
	#"""
	##----------------------------------------------------------------
	#if (CheckEstablishedSibling(clust2_child1_taxa_list, clust2_child2_taxa_list) == False):
		
		## trying the configuration (C, (A, B)) or (B, (A, C))
		#res_clust2_child1_high_leaf_high = CheckR1RelationLevelBased(clust2_child1_High_Level_Taxa_List[0], leaf_taxon)
		#res_clust2_child2_high_leaf_high = CheckR1RelationLevelBased(clust2_child2_High_Level_Taxa_List[0], leaf_taxon)
		
		#if (DEBUG_LEVEL >= 2):
			#fp = open(Output_Text_File, 'a')
			#fp.write('\n ===>>> Obtaining pairwise results ---- ') 
			#fp.write('\n ---- res_clust2_child1_high_leaf_high: ' + str(res_clust2_child1_high_leaf_high))
			#fp.write('\n ---- res_clust2_child2_high_leaf_high: ' + str(res_clust2_child2_high_leaf_high))
			#fp.close()
			
		#"""
		#case 1 - case (C, (A, B))
		#here level based comparison yields (C,A) as 1, and (B,A) as -1
		#"""
		#if (res_clust2_child2_high_leaf_high == 1) and (res_clust2_child1_high_leaf_high != 1):	# == -1):	#sourya
			#src_subtree_node = leaf_mrca_node
			#dest_subtree_node = Curr_tree.mrca(taxon_labels=clust2_child1_taxa_list)
			#Curr_tree = InsertSubTree(Curr_tree, src_subtree_node, dest_subtree_node, Output_Text_File)
			#return Curr_tree
		#"""
		#case 2 - case (B, (A, C))
		#here level based comparison yields (C,A) as -1, and (B,A) as 1
		#"""
		#if (res_clust2_child1_high_leaf_high == 1) and (res_clust2_child2_high_leaf_high != 1):	#== -1):	#sourya
			#src_subtree_node = leaf_mrca_node
			#dest_subtree_node = Curr_tree.mrca(taxon_labels=clust2_child2_taxa_list)
			#Curr_tree = InsertSubTree(Curr_tree, src_subtree_node, dest_subtree_node, Output_Text_File)
			#return Curr_tree
		#"""
		#case 3 - default - (A,(B,C))
		#when none of the above cases are valid
		#"""
		#Curr_tree = MergeSubtrees(Curr_tree, leaf_mrca_node, non_leaf_mrca_node, all_taxa_mrca_node, taxa_list, Output_Text_File)
		#return Curr_tree
	
##-------------------------------------------
#"""
#this function has an input node and it traverses its underlying taxon and 
#lists their level information
#based on that, underlying taxon set is classified
#"""
#def FindLevelTaxa(Curr_tree, src_node, complete_taxa_list):
	#taxa_list = []
	#Level_List = []
	#for n in src_node.leaf_nodes():
		#if n.taxon.label in complete_taxa_list:
			#taxa_list.append(n.taxon.label)
			#Level_List.append(n.level())
			
	#max_level = max(Level_List)
	#min_level = min(Level_List)
	#High_Level_Taxa_List = [taxa_list[i] for i, x in enumerate(Level_List) if x == min_level]
	#Low_Level_Taxa_List = [taxa_list[i] for i, x in enumerate(Level_List) if x == max_level]

	#return taxa_list, High_Level_Taxa_List, Low_Level_Taxa_List
	
##------------------------------------------------------
#"""
#computing average XL information
#"""
#def ComputeAvgXL(inp_taxon, TaxaList, DIST_MAT_TYPE):
	#no_elem = 0
	#sum_xl = 0
	#for x in TaxaList:
		#key1 = (inp_taxon, x)
		#key2 = (x, inp_taxon)
		#if key1 in TaxaPair_Reln_Dict:
			#no_elem = no_elem + 1
			#sum_xl = sum_xl + TaxaPair_Reln_Dict[key1]._GetNormalizedXLSumGeneTrees(DIST_MAT_TYPE)
		#elif key2 in TaxaPair_Reln_Dict:
			#no_elem = no_elem + 1
			#sum_xl = sum_xl + TaxaPair_Reln_Dict[key2]._GetNormalizedXLSumGeneTrees(DIST_MAT_TYPE)

	#if (no_elem == 0):
		#return 0

	#return (sum_xl * 1.0) / no_elem
	
##------------------------------------------------------
#"""
#this function finds the target node onto which one subtree will be merged
#"""
#def FindDestNode(Curr_tree, DestTaxaList, SrcTaxaList, DIST_MAT_TYPE, low_positive):
	#if (len(DestTaxaList) == 1):
		#"""
		#find the sibling internal node and return
		#"""
		#n = Curr_tree.mrca(taxon_labels=DestTaxaList)
		#if (n.parent_node is not None):
			## comment - sourya
			##for x in n.parent_node.child_nodes():
				##if (Node_Label(x) != Node_Label(n)):
					##return x
			## end comment - sourya
			## add - sourya
			#return (n.sister_nodes())[0]
			## end add - sourya
		#else:
			#return n
	#else:
		#Target_Taxa_List = []
		#min_avg_xl = FindAvgXL(DestTaxaList[0], SrcTaxaList, DIST_MAT_TYPE, False)
		##print 'Average XL with ', DestTaxaList[0], '  is ', min_avg_xl
		#min_idx = 0
		#for i in range(1, len(DestTaxaList)):
			#xl = FindAvgXL(DestTaxaList[i], SrcTaxaList, DIST_MAT_TYPE, False) 
			##print 'Average XL with ', DestTaxaList[i], '  is ', xl
			#if (xl < min_avg_xl):
				#min_avg_xl = xl
				#min_idx = i
		#Target_Taxa_List.append(DestTaxaList[min_idx])
		## find the node corresponding to this taxon
		#n = (Curr_tree.mrca(taxon_labels=Target_Taxa_List))
		
		##print 'DestTaxaList: ', DestTaxaList, ' SrcTaxaList: ', SrcTaxaList
		##print 'target node label: ', Node_Label(n)
		
		## two cases - one in which n has a sister leaf node
		## in other case, n has no sister leaf node
		#if (n.parent_node is not None):
			#for x in n.sister_nodes():
				#if (x.is_leaf() == True):
					##print 'Sister node leaf: ', Node_Label(x)
					## a sister leaf node exists
					#return n.parent_node
				#else:
					#for y in x.leaf_nodes(): 
						#if y in DestTaxaList:
							##print 'Sister internal node leaf within DestTaxaList: ', Node_Label(y)
							## a sister internal node has underlying taxon 
							## which is a high level taxon
							#return n.parent_node
			#if (low_positive == True):
				##print 'low positive true: ', Node_Label(n)
				#return n.parent_node
			## otherwise return the sister node
			#return (n.sister_nodes())[0]
		#else:
			#return n

#-------------------------------------------
"""
this function merges both non leaf nodes (taxa cluster) to refine the output species tree
"""
def Merge_Both_NonLeaf(Curr_tree, clust_species_list, idx1, idx2, taxa_list, Output_Text_File, DIST_MAT_TYPE, DIST_MAT_UPDATE):
	
	#----------------------------------------------------------------
	#representing first cluster as (A,B) and second cluster as (C,D)
	#----------------------------------------------------------------
	
	# this is the MRCA node corresponding to the first cluster taxa list
	clust1_mrca_node = Curr_tree.mrca(taxon_labels=clust_species_list[idx1])
	
	# this is the MRCA node corresponding to the second cluster taxa list
	clust2_mrca_node = Curr_tree.mrca(taxon_labels=clust_species_list[idx2])
	
	# MRCA node of the complete set of taxa (it can be original root of multifurcation as well)
	all_taxa_mrca_node = Curr_tree.mrca(taxon_labels=taxa_list)
	
	# child nodes of both these MRCA nodes of respective taxa clusters
	clust1_children = clust1_mrca_node.child_nodes()
	clust1_child1_taxa_list = []
	clust1_child2_taxa_list = []
	# comment - sourya
	#for n in clust1_children[0].leaf_nodes():
		#if n.taxon.label in taxa_list:
			#clust1_child1_taxa_list.append(n.taxon.label)
	# add - sourya
	for n in clust1_children[0].preorder_iter():
		if (n.is_leaf() == True):
			if n.taxon.label in taxa_list:
				clust1_child1_taxa_list.append(n.taxon.label)
	# comment - sourya
	#for n in clust1_children[1].leaf_nodes():
		#if n.taxon.label in taxa_list:
			#clust1_child2_taxa_list.append(n.taxon.label)
	# add - sourya
	for n in clust1_children[1].preorder_iter():
		if (n.is_leaf() == True):
			if n.taxon.label in taxa_list:
				clust1_child2_taxa_list.append(n.taxon.label)

	if (DEBUG_LEVEL >= 2):
		fp = open(Output_Text_File, 'a')
		fp.write('\n ===>>> Merging both non leaf nodes') 
		fp.write('\n ---- First non leaf idx list: ' + str(clust_species_list[idx1]))
		fp.write('\n ---- clust1_child1_taxa_list: ' + str(clust1_child1_taxa_list)) 
		fp.write('\n ---- clust1_child2_taxa_list: ' + str(clust1_child2_taxa_list))
		fp.close()
	
	clust2_children = clust2_mrca_node.child_nodes()
	clust2_child1_taxa_list = []
	clust2_child2_taxa_list = []
	# comment - sourya
	#for n in clust2_children[0].leaf_nodes():
		#if n.taxon.label in taxa_list:
			#clust2_child1_taxa_list.append(n.taxon.label)
	# add - sourya
	for n in clust2_children[0].preorder_iter():
		if (n.is_leaf() == True):
			if n.taxon.label in taxa_list:
				clust2_child1_taxa_list.append(n.taxon.label)
			
	# comment - sourya
	#for n in clust2_children[1].leaf_nodes():
		#if n.taxon.label in taxa_list:
			#clust2_child2_taxa_list.append(n.taxon.label)
	# add - sourya
	for n in clust2_children[1].preorder_iter():
		if (n.is_leaf() == True):
			if n.taxon.label in taxa_list:
				clust2_child2_taxa_list.append(n.taxon.label)

	if (DEBUG_LEVEL >= 2):
		fp = open(Output_Text_File, 'a')
		fp.write('\n ===>>> Merging both non leaf nodes') 
		fp.write('\n ---- Second non leaf idx list: ' + str(clust_species_list[idx2]))
		fp.write('\n ---- clust2_child1_taxa_list: ' + str(clust2_child1_taxa_list)) 
		fp.write('\n ---- clust2_child2_taxa_list: ' + str(clust2_child2_taxa_list))
		fp.close()
	
	#------------------------------
	# obtain XL values for different taxa sets
	# XL(A,B)
	XL_clust1_child1_clust1_child2 = FindAvgXL(clust1_child1_taxa_list, clust1_child2_taxa_list, DIST_MAT_TYPE, True, (DIST_MAT_UPDATE - 1))
	# XL(C,D)
	XL_clust2_child1_clust2_child2 = FindAvgXL(clust2_child1_taxa_list, clust2_child2_taxa_list, DIST_MAT_TYPE, True, (DIST_MAT_UPDATE - 1))
	# XL(A,C)
	XL_clust1_child1_clust2_child1 = FindAvgXL(clust1_child1_taxa_list, clust2_child1_taxa_list, DIST_MAT_TYPE, True, (DIST_MAT_UPDATE - 1))
	# XL(A,D)
	XL_clust1_child1_clust2_child2 = FindAvgXL(clust1_child1_taxa_list, clust2_child2_taxa_list, DIST_MAT_TYPE, True, (DIST_MAT_UPDATE - 1))
	# XL(B,C)
	XL_clust1_child2_clust2_child1 = FindAvgXL(clust1_child2_taxa_list, clust2_child1_taxa_list, DIST_MAT_TYPE, True, (DIST_MAT_UPDATE - 1))
	# XL(B,D)
	XL_clust1_child2_clust2_child2 = FindAvgXL(clust1_child2_taxa_list, clust2_child2_taxa_list, DIST_MAT_TYPE, True, (DIST_MAT_UPDATE - 1))
	
	if (DEBUG_LEVEL >= 2):
		fp = open(Output_Text_File, 'a')
		fp.write('\n ---- XL_clust1_child1_clust1_child2: ' + str(XL_clust1_child1_clust1_child2))
		fp.write('\n ---- XL_clust2_child1_clust2_child2: ' + str(XL_clust2_child1_clust2_child2))
		fp.write('\n ---- XL_clust1_child1_clust2_child1: ' + str(XL_clust1_child1_clust2_child1))
		fp.write('\n ---- XL_clust1_child1_clust2_child2: ' + str(XL_clust1_child1_clust2_child2))
		fp.write('\n ---- XL_clust1_child2_clust2_child1: ' + str(XL_clust1_child2_clust2_child1))
		fp.write('\n ---- XL_clust1_child2_clust2_child2: ' + str(XL_clust1_child2_clust2_child2))
		fp.close()
	
	if (XL_clust1_child1_clust1_child2 <= XL_clust2_child1_clust2_child2):
		"""
		cluster (A,B) will remain intact
		we have to inspect between three configurations
		1) ((A,B),(C,D)), 2) (C, (D, (A, B))), 3) (D, (C, (A, B)))
		"""
		#  XL(C,:) / 2
		XL_clust2_child1_clust1_avg = (XL_clust1_child1_clust2_child1 + XL_clust1_child2_clust2_child1) / 2.0
		# XL(D,:) / 2
		XL_clust2_child2_clust1_avg = (XL_clust1_child1_clust2_child2 + XL_clust1_child2_clust2_child2) / 2.0
		
		"""
		# case 1A - XL(C,D) is lower than both XL(C,:) / 2 and XL(D,:) / 2 (: denotes A and B)
		"""
		if (XL_clust2_child1_clust2_child2 <= XL_clust2_child1_clust1_avg) and (XL_clust2_child1_clust2_child2 <= XL_clust2_child2_clust1_avg):
			# ((A,B),(C,D)) configuration can be employed
			if (DEBUG_LEVEL >= 2):
				fp = open(Output_Text_File, 'a')
				fp.write('\n *** Default Merging Condition ((A,B),(C,D))')
				fp.close()
			Curr_tree = MergeSubtrees(Curr_tree, clust1_mrca_node, clust2_mrca_node, all_taxa_mrca_node, taxa_list, Output_Text_File)
			return Curr_tree
		
		"""
		# case 2A - XL(C,:) / 2 is lower than both XL(C,D) and XL(D,:) / 2 (: denotes A and B)
		"""
		if (XL_clust2_child1_clust1_avg <= XL_clust2_child1_clust2_child2) and (XL_clust2_child1_clust1_avg <= XL_clust2_child2_clust1_avg):
			# configuration (D, (C, (A, B))) can be employed here
			if (DEBUG_LEVEL >= 2):
				fp = open(Output_Text_File, 'a')
				fp.write('\n *** Merging Condition (D, (C, (A, B)))')
				fp.close()
			src_subtree_node = clust1_mrca_node
			dest_subtree_node = Curr_tree.mrca(taxon_labels=clust2_child1_taxa_list)
			Curr_tree = InsertSubTree(Curr_tree, src_subtree_node, dest_subtree_node, Output_Text_File)
			return Curr_tree
		
		"""
		# case 3A - XL(D,:) / 2 is lower than both XL(C,D) and XL(C,:) / 2 (: denotes A and B)
		"""
		if (XL_clust2_child2_clust1_avg <= XL_clust2_child1_clust1_avg) and (XL_clust2_child2_clust1_avg <= XL_clust2_child1_clust2_child2):
			# configuration (C, (D, (A, B))) can be employed here
			if (DEBUG_LEVEL >= 2):
				fp = open(Output_Text_File, 'a')
				fp.write('\n *** Merging Condition (C, (D, (A, B)))')
				fp.close()
			src_subtree_node = clust1_mrca_node
			dest_subtree_node = Curr_tree.mrca(taxon_labels=clust2_child2_taxa_list)
			Curr_tree = InsertSubTree(Curr_tree, src_subtree_node, dest_subtree_node, Output_Text_File)
			return Curr_tree
		
	else:
		"""
		cluster (C,D) will remain intact
		we have to inspect between three configurations
		1) ((A,B),(C,D)), 2) (A, (B, (C, D))), 3) (B, (A, (C, D)))
		"""
		#  XL(A,:) / 2
		XL_clust1_child1_clust2_avg = (XL_clust1_child1_clust2_child1 + XL_clust1_child1_clust2_child2) / 2.0
		# XL(B,:) / 2
		XL_clust1_child2_clust2_avg = (XL_clust1_child2_clust2_child1 + XL_clust1_child2_clust2_child2) / 2.0
		"""
		# case 1B - XL(A,B) is lower than both XL(A,:) / 2 and XL(B,:) / 2 (: denotes C and D)
		"""
		if (XL_clust1_child1_clust1_child2 <= XL_clust1_child1_clust2_avg) and (XL_clust1_child1_clust1_child2 <= XL_clust1_child2_clust2_avg):
			if (DEBUG_LEVEL >= 2):
				fp = open(Output_Text_File, 'a')
				fp.write('\n *** Default Merging Condition ((A,B),(C,D))')
				fp.close()
			Curr_tree = MergeSubtrees(Curr_tree, clust1_mrca_node, clust2_mrca_node, all_taxa_mrca_node, taxa_list, Output_Text_File)
			return Curr_tree
		
		"""
		# case 2B - XL(A,:) / 2 is lower than both XL(A,B) and XL(B,:) / 2 (: denotes C and D)
		"""
		if (XL_clust1_child1_clust2_avg <= XL_clust1_child1_clust1_child2) and (XL_clust1_child1_clust2_avg <= XL_clust1_child2_clust2_avg):
			# the configuration (B, (A, (C, D))) can be employed here
			if (DEBUG_LEVEL >= 2):
				fp = open(Output_Text_File, 'a')
				fp.write('\n *** Merging Condition (B, (A, (C, D)))')
				fp.close()
			src_subtree_node = clust2_mrca_node
			dest_subtree_node = Curr_tree.mrca(taxon_labels=clust1_child1_taxa_list)
			Curr_tree = InsertSubTree(Curr_tree, src_subtree_node, dest_subtree_node, Output_Text_File)
			return Curr_tree
			
		"""
		# case 3B - XL(B,:) / 2 is lower than both XL(A,B) and XL(A,:) / 2 (: denotes C and D)
		"""
		if (XL_clust1_child2_clust2_avg <= XL_clust1_child1_clust1_child2) and (XL_clust1_child2_clust2_avg <= XL_clust1_child1_clust2_avg):
			# the configuration (A, (B, (C, D))) can be employed here
			if (DEBUG_LEVEL >= 2):
				fp = open(Output_Text_File, 'a')
				fp.write('\n *** Merging Condition (A, (B, (C, D)))')
				fp.close()
			src_subtree_node = clust2_mrca_node
			dest_subtree_node = Curr_tree.mrca(taxon_labels=clust1_child2_taxa_list)
			Curr_tree = InsertSubTree(Curr_tree, src_subtree_node, dest_subtree_node, Output_Text_File)
			return Curr_tree
			
	#-----------------------------------------------------
	# comment - sourya

	#"""
	#first we inspect the case when ((A,B),(C,D)) configuration can be employed
	#"""
	#XL_list = [XL_clust1_child1_clust1_child2, XL_clust2_child1_clust2_child2, XL_clust1_child1_clust2_child1, \
		#XL_clust1_child1_clust2_child2, XL_clust1_child2_clust2_child1, XL_clust1_child2_clust2_child2]
	#XL_list.sort()
	
	#"""
	#if the XL count of (A,B) and (C,D) lie within first three indices
	#then merge this cluster pair
	#"""
	#if (XL_list.index(XL_clust1_child1_clust1_child2) <= 2) and (XL_list.index(XL_clust2_child1_clust2_child2) <= 2):
		#if (DEBUG_LEVEL >= 2):
			#fp = open(Output_Text_File, 'a')
			#fp.write('\n *** Default Merging Condition ((A,B),(C,D))')
			#fp.close()
		
		## ((A,B),(C,D)) configuration can be employed
		#Curr_tree = MergeSubtrees(Curr_tree, clust1_mrca_node, clust2_mrca_node, all_taxa_mrca_node, taxa_list, Output_Text_File)
		#return Curr_tree
	##------------------------------------------------------
	#"""
	#otherwise, we have to inspect other configurations
	#"""
	## computing the sum of XL values for all different clusters
	## XL sum for cluster A
	#sum_XL_clust1_child1 = XL_clust1_child1_clust1_child2 + XL_clust1_child1_clust2_child1 + XL_clust1_child1_clust2_child2
	## XL sum for cluster B
	#sum_XL_clust1_child2 = XL_clust1_child1_clust1_child2 + XL_clust1_child2_clust2_child1 + XL_clust1_child2_clust2_child2
	## XL sum for cluster C
	#sum_XL_clust2_child1 = XL_clust2_child1_clust2_child2 + XL_clust1_child1_clust2_child1 + XL_clust1_child2_clust2_child1
	## XL sum for cluster D
	#sum_XL_clust2_child2 = XL_clust2_child1_clust2_child2 + XL_clust1_child1_clust2_child2 + XL_clust1_child2_clust2_child2
	
	#sum_XL_list = [sum_XL_clust1_child1, sum_XL_clust1_child2, sum_XL_clust2_child1, sum_XL_clust2_child2]
		
	#if (DEBUG_LEVEL >= 2):
		#fp = open(Output_Text_File, 'a')
		#fp.write('\n sum_XL_clust1_child1: ' + str(sum_XL_clust1_child1))
		#fp.write('\n sum_XL_clust1_child2: ' + str(sum_XL_clust1_child2))
		#fp.write('\n sum_XL_clust2_child1: ' + str(sum_XL_clust2_child1))
		#fp.write('\n sum_XL_clust2_child2: ' + str(sum_XL_clust2_child2))
		#fp.close()
	
	#if (sum_XL_clust1_child1 == max(sum_XL_list)):
		## the configuration (A, (B, (C, D))) can be employed here
		#if (DEBUG_LEVEL >= 2):
			#fp = open(Output_Text_File, 'a')
			#fp.write('\n *** Merging Condition (A, (B, (C, D)))')
			#fp.close()
		#src_subtree_node = clust2_mrca_node
		#dest_subtree_node = Curr_tree.mrca(taxon_labels=clust1_child2_taxa_list)
		#Curr_tree = InsertSubTree(Curr_tree, src_subtree_node, dest_subtree_node, Output_Text_File)
	#elif (sum_XL_clust1_child2 == max(sum_XL_list)):
		## the configuration (B, (A, (C, D))) can be employed here
		#if (DEBUG_LEVEL >= 2):
			#fp = open(Output_Text_File, 'a')
			#fp.write('\n *** Merging Condition (B, (A, (C, D)))')
			#fp.close()
		#src_subtree_node = clust2_mrca_node
		#dest_subtree_node = Curr_tree.mrca(taxon_labels=clust1_child1_taxa_list)
		#Curr_tree = InsertSubTree(Curr_tree, src_subtree_node, dest_subtree_node, Output_Text_File)
	#elif (sum_XL_clust2_child1 == max(sum_XL_list)):
		## configuration (C, (D, (A, B))) can be employed here
		#if (DEBUG_LEVEL >= 2):
			#fp = open(Output_Text_File, 'a')
			#fp.write('\n *** Merging Condition (C, (D, (A, B)))')
			#fp.close()
		#src_subtree_node = clust1_mrca_node
		#dest_subtree_node = Curr_tree.mrca(taxon_labels=clust2_child2_taxa_list)
		#Curr_tree = InsertSubTree(Curr_tree, src_subtree_node, dest_subtree_node, Output_Text_File)
	#else:	#if (sum_XL_clust2_child2 == max(sum_XL_list)):
		## configuration (D, (C, (A, B))) can be employed here
		#if (DEBUG_LEVEL >= 2):
			#fp = open(Output_Text_File, 'a')
			#fp.write('\n *** Merging Condition (D, (C, (A, B)))')
			#fp.close()
		#src_subtree_node = clust1_mrca_node
		#dest_subtree_node = Curr_tree.mrca(taxon_labels=clust2_child1_taxa_list)
		#Curr_tree = InsertSubTree(Curr_tree, src_subtree_node, dest_subtree_node, Output_Text_File)
		
	#return Curr_tree
	
	# end comment - sourya
	#-----------------------------------------------------

##-------------------------------------------
#"""
#this function merges both non leaf nodes (taxa cluster) to refine the output species tree
#"""
#def Merge_Both_NonLeaf_OLD(Curr_tree, clust_species_list, idx1, idx2, taxa_list, Output_Text_File, DIST_MAT_TYPE):
	
	## this is the MRCA node corresponding to the first cluster taxa list
	#clust1_mrca_node = Curr_tree.mrca(taxon_labels=clust_species_list[idx1])
	
	## this is the MRCA node corresponding to the second cluster taxa list
	#clust2_mrca_node = Curr_tree.mrca(taxon_labels=clust_species_list[idx2])
	
	## MRCA node of the complete set of taxa (it can be original root of multifurcation as well)
	#all_taxa_mrca_node = Curr_tree.mrca(taxon_labels=taxa_list)
	
	## child nodes of both these MRCA nodes of respective taxa clusters
	#clust1_children = clust1_mrca_node.child_nodes()
	#clust1_child1_taxa_list, clust1_child1_High_Level_Taxa_List, clust1_child1_Low_Level_Taxa_List = \
		#FindLevelTaxa(Curr_tree, clust1_children[0], taxa_list)
	#clust1_child2_taxa_list, clust1_child2_High_Level_Taxa_List, clust1_child2_Low_Level_Taxa_List = \
		#FindLevelTaxa(Curr_tree, clust1_children[1], taxa_list)
	
	#if (DEBUG_LEVEL >= 2):
		#fp = open(Output_Text_File, 'a')
		#fp.write('\n ===>>> Merging both non leaf nodes') 
		#fp.write('\n ---- First non leaf idx list: ' + str(clust_species_list[idx1]))
		#fp.write('\n ---- clust1_child1_High_Level_Taxa_List: ' + str(clust1_child1_High_Level_Taxa_List)) 
		#fp.write('\n ---- clust1_child1_Low_Level_Taxa_List: ' + str(clust1_child1_Low_Level_Taxa_List))
		#fp.write('\n ---- clust1_child2_High_Level_Taxa_List: ' + str(clust1_child2_High_Level_Taxa_List)) 
		#fp.write('\n ---- clust1_child2_Low_Level_Taxa_List: ' + str(clust1_child2_Low_Level_Taxa_List))
		#fp.close()
	
	#clust2_children = clust2_mrca_node.child_nodes()
	#clust2_child1_taxa_list, clust2_child1_High_Level_Taxa_List, clust2_child1_Low_Level_Taxa_List = \
		#FindLevelTaxa(Curr_tree, clust2_children[0], taxa_list)
	#clust2_child2_taxa_list, clust2_child2_High_Level_Taxa_List, clust2_child2_Low_Level_Taxa_List = \
		#FindLevelTaxa(Curr_tree, clust2_children[1], taxa_list)
	
	#if (DEBUG_LEVEL >= 2):
		#fp = open(Output_Text_File, 'a')
		#fp.write('\n ===>>> Merging both non leaf nodes') 
		#fp.write('\n ---- Second non leaf idx list: ' + str(clust_species_list[idx2]))
		#fp.write('\n ---- clust2_child1_High_Level_Taxa_List: ' + str(clust2_child1_High_Level_Taxa_List)) 
		#fp.write('\n ---- clust2_child1_Low_Level_Taxa_List: ' + str(clust2_child1_Low_Level_Taxa_List))
		#fp.write('\n ---- clust2_child2_High_Level_Taxa_List: ' + str(clust2_child2_High_Level_Taxa_List)) 
		#fp.write('\n ---- clust2_child2_Low_Level_Taxa_List: ' + str(clust2_child2_Low_Level_Taxa_List))
		#fp.close()
	
	##----------------------------------------------------------------
	#""" 
	#representing first cluster as (A,B) and second cluster as (C,D)
	#following operations can be performed only if A and B are not established sibling couplet
	#"""
	##----------------------------------------------------------------
	#if (CheckEstablishedSibling(clust1_child1_taxa_list, clust1_child2_taxa_list) == False):
		## add - sourya
		
		#res_clust1_child1_high_clust2_child1_high = \
			#Check_R1Reln_Majority(clust1_child1_High_Level_Taxa_List[0], clust2_child1_High_Level_Taxa_List[0])

		#res_clust1_child1_high_clust2_child2_high = \
			#Check_R1Reln_Majority(clust1_child1_High_Level_Taxa_List[0], clust2_child2_High_Level_Taxa_List[0])

		#res_clust1_child2_high_clust2_child1_high = \
			#Check_R1Reln_Majority(clust1_child2_High_Level_Taxa_List[0], clust2_child1_High_Level_Taxa_List[0])

		#res_clust1_child2_high_clust2_child2_high = \
			#Check_R1Reln_Majority(clust1_child2_High_Level_Taxa_List[0], clust2_child2_High_Level_Taxa_List[0])
		
		
		## debug - sourya
		#res_clust1_child1_low_clust2_child1_high = \
			#Check_R1Reln_Majority(clust1_child1_Low_Level_Taxa_List[0], clust2_child1_High_Level_Taxa_List[0])

		#res_clust1_child1_low_clust2_child2_high = \
			#Check_R1Reln_Majority(clust1_child1_Low_Level_Taxa_List[0], clust2_child2_High_Level_Taxa_List[0])

		#res_clust1_child2_low_clust2_child1_high = \
			#Check_R1Reln_Majority(clust1_child2_Low_Level_Taxa_List[0], clust2_child1_High_Level_Taxa_List[0])

		#res_clust1_child2_low_clust2_child2_high = \
			#Check_R1Reln_Majority(clust1_child2_Low_Level_Taxa_List[0], clust2_child2_High_Level_Taxa_List[0])
		
		#if (DEBUG_LEVEL >= 2):
			#fp = open(Output_Text_File, 'a')
			#fp.write('\n ===>>> Obtaining pairwise results ---- ') 
			#fp.write('\n ---- res_clust1_child1_high_clust2_child1_high: ' + str(res_clust1_child1_high_clust2_child1_high))
			#fp.write('\n ---- res_clust1_child1_high_clust2_child2_high: ' + str(res_clust1_child1_high_clust2_child2_high))
			#fp.write('\n ---- res_clust1_child2_high_clust2_child1_high: ' + str(res_clust1_child2_high_clust2_child1_high))
			#fp.write('\n ---- res_clust1_child2_high_clust2_child2_high: ' + str(res_clust1_child2_high_clust2_child2_high))
			#fp.write('\n ---- res_clust1_child1_low_clust2_child1_high: ' + str(res_clust1_child1_low_clust2_child1_high))
			#fp.write('\n ---- res_clust1_child1_low_clust2_child2_high: ' + str(res_clust1_child1_low_clust2_child2_high))
			#fp.write('\n ---- res_clust1_child2_low_clust2_child1_high: ' + str(res_clust1_child2_low_clust2_child1_high))
			#fp.write('\n ---- res_clust1_child2_low_clust2_child2_high: ' + str(res_clust1_child2_low_clust2_child2_high))
			#fp.close()
			
		##------------------------------------------
		#res_clust1_child1_low_clust2_child1_low = \
			#Check_R1Reln_Majority(clust1_child1_Low_Level_Taxa_List[0], clust2_child1_Low_Level_Taxa_List[0])

		#res_clust1_child1_low_clust2_child2_low = \
			#Check_R1Reln_Majority(clust1_child1_Low_Level_Taxa_List[0], clust2_child2_Low_Level_Taxa_List[0])

		#res_clust1_child2_low_clust2_child1_low = \
			#Check_R1Reln_Majority(clust1_child2_Low_Level_Taxa_List[0], clust2_child1_Low_Level_Taxa_List[0])

		#res_clust1_child2_low_clust2_child2_low = \
			#Check_R1Reln_Majority(clust1_child2_Low_Level_Taxa_List[0], clust2_child2_Low_Level_Taxa_List[0])

		#if (DEBUG_LEVEL >= 2):
			#fp = open(Output_Text_File, 'a')
			#fp.write('\n ===>>> Obtaining pairwise results ---- ') 
			#fp.write('\n ---- res_clust1_child1_low_clust2_child1_low: ' + str(res_clust1_child1_low_clust2_child1_low))
			#fp.write('\n ---- res_clust1_child1_low_clust2_child2_low: ' + str(res_clust1_child1_low_clust2_child2_low))
			#fp.write('\n ---- res_clust1_child2_low_clust2_child1_low: ' + str(res_clust1_child2_low_clust2_child1_low))
			#fp.write('\n ---- res_clust1_child2_low_clust2_child2_low: ' + str(res_clust1_child2_low_clust2_child2_low))
			#fp.close()

		#if (res_clust1_child1_low_clust2_child1_low >= 1) and (res_clust1_child1_low_clust2_child2_low >= 1) and \
			#(res_clust1_child2_low_clust2_child1_low >= 1) and (res_clust1_child2_low_clust2_child2_low >= 1):
			#low_positive = True
		#else:
			#low_positive = False
		##------------------------------------------
		#if (res_clust1_child1_high_clust2_child1_high == 1) and (res_clust1_child1_high_clust2_child2_high == 1) and \
			#(((res_clust1_child2_high_clust2_child1_high == 2) and (res_clust1_child2_high_clust2_child2_high == 2)) or \
				#((res_clust1_child2_high_clust2_child1_high == 1) and (res_clust1_child2_high_clust2_child2_high == 2)) or 
			#((res_clust1_child2_high_clust2_child1_high == 2) and (res_clust1_child2_high_clust2_child2_high == 1))):
			
			#if (DEBUG_LEVEL >= 2):
				#fp = open(Output_Text_File, 'a')
				#fp.write('\n *** Condition 1A')
				#fp.close()
			#src_subtree_node = clust2_mrca_node
			
			## this is the taxa list which will be re-positioned
			#Complete_Src_Taxa_List = []
			#for x in clust2_child1_High_Level_Taxa_List:
				#Complete_Src_Taxa_List.append(x)
			#for x in clust2_child2_High_Level_Taxa_List:
				#Complete_Src_Taxa_List.append(x)
				
			## this is the taxa list which will be checked for repositioning the earlier mentioned taxa list
			#Complete_Dest_Taxa_List = []
			#for x in clust1_child1_High_Level_Taxa_List:
				#Complete_Dest_Taxa_List.append(x)
			#for x in clust1_child2_High_Level_Taxa_List:
				#Complete_Dest_Taxa_List.append(x)
			
			#dest_subtree_node = FindDestNode(Curr_tree, Complete_Dest_Taxa_List, Complete_Src_Taxa_List, DIST_MAT_TYPE, low_positive)
			
			#Curr_tree = InsertSubTree(Curr_tree, src_subtree_node, dest_subtree_node, Output_Text_File)
			#return Curr_tree

		#if (((res_clust1_child1_high_clust2_child1_high == 2) and (res_clust1_child1_high_clust2_child2_high == 2)) or \
			#((res_clust1_child1_high_clust2_child1_high == 1) and (res_clust1_child1_high_clust2_child2_high == 2)) or \
				#((res_clust1_child1_high_clust2_child1_high == 2) and (res_clust1_child1_high_clust2_child2_high == 1))) and \
			#(res_clust1_child2_high_clust2_child1_high == 1) and (res_clust1_child2_high_clust2_child2_high == 1):
			
			#if (DEBUG_LEVEL >= 2):
				#fp = open(Output_Text_File, 'a')
				#fp.write('\n *** Condition 2A')
				#fp.close()
			#src_subtree_node = clust2_mrca_node
			
			## this is the taxa list which will be re-positioned
			#Complete_Src_Taxa_List = []
			#for x in clust2_child1_High_Level_Taxa_List:
				#Complete_Src_Taxa_List.append(x)
			#for x in clust2_child2_High_Level_Taxa_List:
				#Complete_Src_Taxa_List.append(x)
				
			## this is the taxa list which will be checked for repositioning the earlier mentioned taxa list
			#Complete_Dest_Taxa_List = []
			#for x in clust1_child1_High_Level_Taxa_List:
				#Complete_Dest_Taxa_List.append(x)
			#for x in clust1_child2_High_Level_Taxa_List:
				#Complete_Dest_Taxa_List.append(x)
			
			#dest_subtree_node = FindDestNode(Curr_tree, Complete_Dest_Taxa_List, Complete_Src_Taxa_List, DIST_MAT_TYPE, low_positive)

			#Curr_tree = InsertSubTree(Curr_tree, src_subtree_node, dest_subtree_node, Output_Text_File)
			#return Curr_tree

		#if (res_clust1_child1_high_clust2_child1_high == 1) and (res_clust1_child1_high_clust2_child2_high == 1) and \
			#(res_clust1_child2_high_clust2_child1_high == 1) and (res_clust1_child2_high_clust2_child2_high == 1):

			## configuration (B, (A, (C, D))) or (A, (B, (C, D))) can be satisfied
			
			## we check if all low values are 1 - in such a case, we proceed with simple merging operation
			## otherwise, we follow earlier operation
			#if (res_clust1_child1_low_clust2_child1_low == 1) and (res_clust1_child1_low_clust2_child2_low == 1) and \
				#(res_clust1_child2_low_clust2_child1_low == 1) and (res_clust1_child2_low_clust2_child2_low == 1):
				
				#if (DEBUG_LEVEL >= 2):
					#fp = open(Output_Text_File, 'a')
					#fp.write('\n *** Condition 3A')
					#fp.close()
				
				## simple merging operation
				## configuration ((A,B),(C,D)) will be employed
				#Curr_tree = MergeSubtrees(Curr_tree, clust1_mrca_node, clust2_mrca_node, all_taxa_mrca_node, taxa_list, Output_Text_File)
				#return Curr_tree
			
			#else:
				
				#if (DEBUG_LEVEL >= 2):
					#fp = open(Output_Text_File, 'a')
					#fp.write('\n *** Condition 4A')
					#fp.close()
				#src_subtree_node = clust2_mrca_node
				
				## this is the taxa list which will be re-positioned
				#Complete_Src_Taxa_List = []
				#for x in clust2_child1_High_Level_Taxa_List:
					#Complete_Src_Taxa_List.append(x)
				#for x in clust2_child2_High_Level_Taxa_List:
					#Complete_Src_Taxa_List.append(x)
					
				## this is the taxa list which will be checked for repositioning the earlier mentioned taxa list
				#Complete_Dest_Taxa_List = []
				#for x in clust1_child1_High_Level_Taxa_List:
					#Complete_Dest_Taxa_List.append(x)
				#for x in clust1_child2_High_Level_Taxa_List:
					#Complete_Dest_Taxa_List.append(x)
				
				#dest_subtree_node = FindDestNode(Curr_tree, Complete_Dest_Taxa_List, Complete_Src_Taxa_List, DIST_MAT_TYPE, low_positive)
				
				#Curr_tree = InsertSubTree(Curr_tree, src_subtree_node, dest_subtree_node, Output_Text_File)
				#return Curr_tree
				
	##----------------------------------------------------------------
	#""" 
	#representing first cluster as (A,B) and second cluster as (C,D)
	#following operations can be performed only if C and D are not established sibling couplet
	#"""
	##----------------------------------------------------------------
	#if (CheckEstablishedSibling(clust2_child1_taxa_list, clust2_child2_taxa_list) == False):
	
		#res_clust2_child1_high_clust1_child1_high = \
			#Check_R1Reln_Majority(clust2_child1_High_Level_Taxa_List[0], clust1_child1_High_Level_Taxa_List[0])

		#res_clust2_child1_high_clust1_child2_high = \
			#Check_R1Reln_Majority(clust2_child1_High_Level_Taxa_List[0], clust1_child2_High_Level_Taxa_List[0])

		#res_clust2_child2_high_clust1_child1_high = \
			#Check_R1Reln_Majority(clust2_child2_High_Level_Taxa_List[0], clust1_child1_High_Level_Taxa_List[0])

		#res_clust2_child2_high_clust1_child2_high = \
			#Check_R1Reln_Majority(clust2_child2_High_Level_Taxa_List[0], clust1_child2_High_Level_Taxa_List[0])


		## debug - sourya
		#res_clust2_child1_low_clust1_child1_high = \
			#Check_R1Reln_Majority(clust2_child1_Low_Level_Taxa_List[0], clust1_child1_High_Level_Taxa_List[0])

		#res_clust2_child1_low_clust1_child2_high = \
			#Check_R1Reln_Majority(clust2_child1_Low_Level_Taxa_List[0], clust1_child2_High_Level_Taxa_List[0])

		#res_clust2_child2_low_clust1_child1_high = \
			#Check_R1Reln_Majority(clust2_child2_Low_Level_Taxa_List[0], clust1_child1_High_Level_Taxa_List[0])

		#res_clust2_child2_low_clust1_child2_high = \
			#Check_R1Reln_Majority(clust2_child2_Low_Level_Taxa_List[0], clust1_child2_High_Level_Taxa_List[0])
		

		#if (DEBUG_LEVEL >= 2):
			#fp = open(Output_Text_File, 'a')
			#fp.write('\n ===>>> Obtaining pairwise results ---- ') 
			#fp.write('\n ---- res_clust2_child1_high_clust1_child1_high: ' + str(res_clust2_child1_high_clust1_child1_high))
			#fp.write('\n ---- res_clust2_child1_high_clust1_child2_high: ' + str(res_clust2_child1_high_clust1_child2_high))
			#fp.write('\n ---- res_clust2_child2_high_clust1_child1_high: ' + str(res_clust2_child2_high_clust1_child1_high))
			#fp.write('\n ---- res_clust2_child2_high_clust1_child2_high: ' + str(res_clust2_child2_high_clust1_child2_high))
			#fp.write('\n ---- res_clust2_child1_low_clust1_child1_high: ' + str(res_clust2_child1_low_clust1_child1_high))
			#fp.write('\n ---- res_clust2_child1_low_clust1_child2_high: ' + str(res_clust2_child1_low_clust1_child2_high))
			#fp.write('\n ---- res_clust2_child2_low_clust1_child1_high: ' + str(res_clust2_child2_low_clust1_child1_high))
			#fp.write('\n ---- res_clust2_child2_low_clust1_child2_high: ' + str(res_clust2_child2_low_clust1_child2_high))
			#fp.close()
			
		##------------------------------------------
		#res_clust2_child1_low_clust1_child1_low = \
			#Check_R1Reln_Majority(clust2_child1_Low_Level_Taxa_List[0], clust1_child1_Low_Level_Taxa_List[0])

		#res_clust2_child1_low_clust1_child2_low = \
			#Check_R1Reln_Majority(clust2_child1_Low_Level_Taxa_List[0], clust1_child2_Low_Level_Taxa_List[0])

		#res_clust2_child2_low_clust1_child1_low = \
			#Check_R1Reln_Majority(clust2_child2_Low_Level_Taxa_List[0], clust1_child1_Low_Level_Taxa_List[0])

		#res_clust2_child2_low_clust1_child2_low = \
			#Check_R1Reln_Majority(clust2_child2_Low_Level_Taxa_List[0], clust1_child2_Low_Level_Taxa_List[0])

		#if (DEBUG_LEVEL >= 2):
			#fp = open(Output_Text_File, 'a')
			#fp.write('\n ===>>> Obtaining pairwise results ---- ') 
			#fp.write('\n ---- res_clust2_child1_low_clust1_child1_low: ' + str(res_clust2_child1_low_clust1_child1_low))
			#fp.write('\n ---- res_clust2_child1_low_clust1_child2_low: ' + str(res_clust2_child1_low_clust1_child2_low))
			#fp.write('\n ---- res_clust2_child2_low_clust1_child1_low: ' + str(res_clust2_child2_low_clust1_child1_low))
			#fp.write('\n ---- res_clust2_child2_low_clust1_child2_low: ' + str(res_clust2_child2_low_clust1_child2_low))
			#fp.close()
			
		#if (res_clust2_child1_low_clust1_child1_low >= 1) and (res_clust2_child1_low_clust1_child2_low >= 1) and \
			#(res_clust2_child2_low_clust1_child1_low >= 1) and (res_clust2_child2_low_clust1_child2_low >= 1):
			#low_positive = True
		#else:
			#low_positive = False

		##------------------------------------------

		#if (res_clust2_child1_high_clust1_child1_high == 1) and (res_clust2_child1_high_clust1_child2_high == 1) and \
			#(((res_clust2_child2_high_clust1_child1_high == 2) and (res_clust2_child2_high_clust1_child2_high == 2)) or \
				#((res_clust2_child2_high_clust1_child1_high == 1) and (res_clust2_child2_high_clust1_child2_high == 2)) or \
					#((res_clust2_child2_high_clust1_child1_high == 2) and (res_clust2_child2_high_clust1_child2_high == 1))):

			#if (DEBUG_LEVEL >= 2):
				#fp = open(Output_Text_File, 'a')
				#fp.write('\n *** Condition 1B')
				#fp.close()
			#src_subtree_node = clust1_mrca_node
			
			## this is the taxa list which will be re-positioned
			#Complete_Src_Taxa_List = []
			#for x in clust1_child1_High_Level_Taxa_List:
				#Complete_Src_Taxa_List.append(x)
			#for x in clust1_child2_High_Level_Taxa_List:
				#Complete_Src_Taxa_List.append(x)
				
			## this is the taxa list which will be checked for repositioning the earlier mentioned taxa list
			#Complete_Dest_Taxa_List = []
			#for x in clust2_child1_High_Level_Taxa_List:
				#Complete_Dest_Taxa_List.append(x)
			#for x in clust2_child2_High_Level_Taxa_List:
				#Complete_Dest_Taxa_List.append(x)
			
			#dest_subtree_node = FindDestNode(Curr_tree, Complete_Dest_Taxa_List, Complete_Src_Taxa_List, DIST_MAT_TYPE, low_positive)
			
			#Curr_tree = InsertSubTree(Curr_tree, src_subtree_node, dest_subtree_node, Output_Text_File)
			#return Curr_tree
	
		#if (((res_clust2_child1_high_clust1_child1_high == 2) and (res_clust2_child1_high_clust1_child2_high == 2)) or \
			#((res_clust2_child1_high_clust1_child1_high == 2) and (res_clust2_child1_high_clust1_child2_high == 1)) or \
				#((res_clust2_child1_high_clust1_child1_high == 1) and (res_clust2_child1_high_clust1_child2_high == 2))) and \
			#(res_clust2_child2_high_clust1_child1_high == 1) and (res_clust2_child2_high_clust1_child2_high == 1):
				
			#if (DEBUG_LEVEL >= 2):
				#fp = open(Output_Text_File, 'a')
				#fp.write('\n *** Condition 2B')
				#fp.close()
			#src_subtree_node = clust1_mrca_node
			
			## this is the taxa list which will be re-positioned
			#Complete_Src_Taxa_List = []
			#for x in clust1_child1_High_Level_Taxa_List:
				#Complete_Src_Taxa_List.append(x)
			#for x in clust1_child2_High_Level_Taxa_List:
				#Complete_Src_Taxa_List.append(x)
				
			## this is the taxa list which will be checked for repositioning the earlier mentioned taxa list
			#Complete_Dest_Taxa_List = []
			#for x in clust2_child1_High_Level_Taxa_List:
				#Complete_Dest_Taxa_List.append(x)
			#for x in clust2_child2_High_Level_Taxa_List:
				#Complete_Dest_Taxa_List.append(x)
			
			#dest_subtree_node = FindDestNode(Curr_tree, Complete_Dest_Taxa_List, Complete_Src_Taxa_List, DIST_MAT_TYPE, low_positive)
			
			#Curr_tree = InsertSubTree(Curr_tree, src_subtree_node, dest_subtree_node, Output_Text_File)
			#return Curr_tree

		
		#if (res_clust2_child1_high_clust1_child1_high == 1) and (res_clust2_child1_high_clust1_child2_high == 1) and \
			#(res_clust2_child2_high_clust1_child1_high == 1) and (res_clust2_child2_high_clust1_child2_high == 1):
			
			## configuration (C, (D, (A, B))) or (D, (C, (A, B)) can be satisfied
			
			## we check if all low values are 1 - in such a case, we proceed with simple merging operation
			## otherwise, we follow earlier operation
			#if (res_clust2_child1_low_clust1_child1_low == 1) and (res_clust2_child1_low_clust1_child2_low == 1) and \
				#(res_clust2_child2_low_clust1_child1_low == 1) and (res_clust2_child2_low_clust1_child2_low == 1):
				
				#if (DEBUG_LEVEL >= 2):
					#fp = open(Output_Text_File, 'a')
					#fp.write('\n *** Condition 3B')
					#fp.close()
				
				## simple merging operation
				## configuration ((A,B),(C,D)) will be employed
				#Curr_tree = MergeSubtrees(Curr_tree, clust1_mrca_node, clust2_mrca_node, all_taxa_mrca_node, taxa_list, Output_Text_File)
				#return Curr_tree
			
			#else:
				
				#if (DEBUG_LEVEL >= 2):
					#fp = open(Output_Text_File, 'a')
					#fp.write('\n *** Condition 4B')
					#fp.close()
				#src_subtree_node = clust1_mrca_node
				
				## this is the taxa list which will be re-positioned
				#Complete_Src_Taxa_List = []
				#for x in clust1_child1_High_Level_Taxa_List:
					#Complete_Src_Taxa_List.append(x)
				#for x in clust1_child2_High_Level_Taxa_List:
					#Complete_Src_Taxa_List.append(x)
					
				## this is the taxa list which will be checked for repositioning the earlier mentioned taxa list
				#Complete_Dest_Taxa_List = []
				#for x in clust2_child1_High_Level_Taxa_List:
					#Complete_Dest_Taxa_List.append(x)
				#for x in clust2_child2_High_Level_Taxa_List:
					#Complete_Dest_Taxa_List.append(x)
				
				#dest_subtree_node = FindDestNode(Curr_tree, Complete_Dest_Taxa_List, Complete_Src_Taxa_List, DIST_MAT_TYPE, low_positive)
				
				#Curr_tree = InsertSubTree(Curr_tree, src_subtree_node, dest_subtree_node, Output_Text_File)
				#return Curr_tree
			
	## otherwise the configuration ((A,B),(C,D)) will be employed
	#Curr_tree = MergeSubtrees(Curr_tree, clust1_mrca_node, clust2_mrca_node, all_taxa_mrca_node, taxa_list, Output_Text_File)
	
	#return Curr_tree	

#-------------------------------------------
"""
this is a new function to merge taxa clusters for agglomerative clustering - sourya
"""
def Merge_Cluster_Pair_New(Curr_tree, clust_species_list, min_idx_i, min_idx_j, taxa_list, Output_Text_File, DIST_MAT_TYPE, DIST_MAT_UPDATE):

	isleaf_clust1 = IsLeafCluster(clust_species_list, min_idx_i)
	isleaf_clust2 = IsLeafCluster(clust_species_list, min_idx_j)
	if (isleaf_clust1 == True) and (isleaf_clust2 == True):
		Curr_tree = Merge_Leaves(Curr_tree, clust_species_list, min_idx_i, min_idx_j, taxa_list, Output_Text_File)
	elif (isleaf_clust1 == True) and (isleaf_clust2 == False):
		Curr_tree = Merge_Leaf_NonLeaf(Curr_tree, clust_species_list, min_idx_i, min_idx_j, taxa_list, Output_Text_File, DIST_MAT_TYPE, DIST_MAT_UPDATE)
	elif (isleaf_clust1 == False) and (isleaf_clust2 == True):
		Curr_tree = Merge_Leaf_NonLeaf(Curr_tree, clust_species_list, min_idx_j, min_idx_i, taxa_list, Output_Text_File, DIST_MAT_TYPE, DIST_MAT_UPDATE)
	elif (isleaf_clust1 == False) and (isleaf_clust2 == False):
		Curr_tree = Merge_Both_NonLeaf(Curr_tree, clust_species_list, min_idx_i, min_idx_j, taxa_list, Output_Text_File, DIST_MAT_TYPE, DIST_MAT_UPDATE)
	
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
def ResolveMultifurcation(Curr_tree, clust_species_list, no_of_input_clusters, Output_Text_File, \
	NJ_RULE_USED, DIST_MAT_TYPE, DIST_MAT_UPDATE, NJ_MERGE_CLUST):
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

	# comment - sourya
	#---------------------------------------
	# using single taxon as a representative of the taxa cluster
	#Fill_DistMat_SingleEntry(DistMat, no_of_clust, clust_species_list, DIST_MAT_TYPE)
	
	# using average information from a taxa cluster
	Fill_DistMat_AvgEntry(DistMat, no_of_clust, clust_species_list, DIST_MAT_TYPE, (DIST_MAT_UPDATE - 1))
	#---------------------------------------
	# end comment - sourya
	
	# loop to execute the agglomerative clustering
	while(no_of_clust > 2):
		# add - sourya
		"""
		instead of approximating the cluster contents
		we use the XL value computed from individual cluster elements
		"""
		#Fill_DistMat_SingleEntry(DistMat, no_of_clust, clust_species_list, DIST_MAT_TYPE, (DIST_MAT_UPDATE - 1))
		#Fill_DistMat_AvgEntry(DistMat, no_of_clust, clust_species_list, DIST_MAT_TYPE, (DIST_MAT_UPDATE - 1))
		# end add - sourya
		
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
		if (NJ_MERGE_CLUST == 1):
			Curr_tree = Merge_Cluster_Pair(Curr_tree, clust_species_list, min_idx_i, min_idx_j, taxa_list, Output_Text_File)
		else:
			Curr_tree = Merge_Cluster_Pair_New(Curr_tree, clust_species_list, min_idx_i, min_idx_j, taxa_list, Output_Text_File, DIST_MAT_TYPE, DIST_MAT_UPDATE)
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
		
		# comment - sourya
		# now fill the elements of the new added row and column
		for k in range(no_of_clust):
			# for any index k, the number of extra lineage is the maximum of these three quantities
			if (NJ_RULE_USED == AGGLO_CLUST):
				if (DIST_MAT_UPDATE == 1):
					DistMat[k][no_of_clust] = (DistMat[k][min_idx_i] + DistMat[k][min_idx_j]) / 2.0
				elif (DIST_MAT_UPDATE == 2):
					DistMat[k][no_of_clust] = min(DistMat[k][min_idx_i], DistMat[k][min_idx_j])	# modified - sourya
				else:
					DistMat[k][no_of_clust] = max(DistMat[k][min_idx_i], DistMat[k][min_idx_j])	# modified - sourya
				#elif (DIST_MAT_UPDATE == 3):
					#DistMat[k][no_of_clust] = (DistMat[k][min_idx_i] * len(clust_species_list[min_idx_i]) + \
						#DistMat[k][min_idx_j] * len(clust_species_list[min_idx_j])) / (len(clust_species_list[min_idx_i]) + len(clust_species_list[min_idx_j]))
			else:
				DistMat[k][no_of_clust] = (DistMat[k][min_idx_i] + DistMat[k][min_idx_j] - DistMat[min_idx_i][min_idx_j]) / 2.0
			# symmetric property
			DistMat[no_of_clust][k] = DistMat[k][no_of_clust]
		# end comment - sourya
		
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
		
		# comment - sourya
		#clust_species_list.append(taxa_list)    
		
		# add - sourya
		taxa_list_mrca_node = Curr_tree.mrca(taxon_labels=taxa_list)
		preorder_taxa_list = []
		for n in taxa_list_mrca_node.preorder_iter():
			if (n.is_leaf() == True):
				if n.taxon.label in taxa_list:
					preorder_taxa_list.append(n.taxon.label)
		clust_species_list.append(preorder_taxa_list)   
		# end add - sourya
		
		# decrement the number of clusters considered
		no_of_clust = no_of_clust - 1
	
	return
            
#-------------------------------------------
# this function refines input supertree such that the supertree becomes binary
# this is required for proper benchmarking with existing binary tree construction methods on 
# ILS sorting
def Refine_Supertree_Binary_Form(Curr_tree, Output_Text_File, NJ_RULE_USED, DIST_MAT_TYPE, DIST_MAT_UPDATE, NJ_MERGE_CLUST):

	# comment - sourya
	#----------------------------------------
	## add - sourya
	## maintain the list of couplets
	#for curr_node in Curr_tree.postorder_internal_node_iter():
		#curr_node_children = curr_node.child_nodes()
		#if (len(curr_node_children) == 2) and (curr_node_children[0].is_leaf() == True) and (curr_node_children[1].is_leaf() == True):
			## add - sourya
			#clust_idx_child1 = Taxa_Info_Dict[curr_node_children[0].taxon.label]._Get_Taxa_Part_Clust_Idx()
			#clust_idx_child2 = Taxa_Info_Dict[curr_node_children[1].taxon.label]._Get_Taxa_Part_Clust_Idx()
			#if (clust_idx_child1 == clust_idx_child2):	# condition add - sourya
				#subl = []
				#for x in curr_node_children:
					#if (x.is_leaf() == True):
						#subl.append(x.taxon.label)
				#Sibling_Couplet_List.append(subl)
	## end add - sourya
	
	#if (DEBUG_LEVEL >= 2):
		#fp = open(Output_Text_File, 'a')
		#fp.write('\n Sibling_Couplet_List ' + str(Sibling_Couplet_List))
		#fp.close()
	#----------------------------------------
	
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
			
			# comment - sourya
			## examine individual nodes of the current node's children list
			#for x in curr_node_children:
				#if (x.is_leaf() == True):
					#subl = []
					#subl.append(x.taxon.label)
				#else:
					#subl = GetTaxaUnderInternalNode(x)
				#clust_species_list.append(subl)
			# end comment - sourya

			# add - sourya
			for x in curr_node_children:
				subl = []
				for n in x.preorder_iter():
					if (n.is_leaf() == True):
						subl.append(n.taxon.label)
				clust_species_list.append(subl)
			# end add - sourya
			
			# call the resolving routine
			ResolveMultifurcation(Curr_tree, clust_species_list, len(curr_node_children), Output_Text_File, \
				NJ_RULE_USED, DIST_MAT_TYPE, DIST_MAT_UPDATE, NJ_MERGE_CLUST)
