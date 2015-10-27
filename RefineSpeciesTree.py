#!/usr/bin/env python

import Header
from Header import *
import UtilFunc
from UtilFunc import *

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

#-------------------------------------------
"""
this function processes one internal node (basically the children list)
to resolve multifurcation
"""
def ResolveMultifurcation(Gene_TreeList, Curr_tree, clust_species_list, no_of_input_clusters, Output_Text_File):
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
	# add - sourya
	for i in range(no_of_clust - 1):
		for j in range(i+1, no_of_clust):
			# here both clust_species_list[i] and clust_species_list[j]
			# are one element lists (according to their construction)
			# we have extracted the corresponding element by using [0] operator (extracting first element)
			x1 = clust_species_list[i][0]
			x2 = clust_species_list[j][0]
			key1 = (x1, x2)
			key2 = (x2, x1)
			#print 'key1: ', key1, ' key2: ', key2
			if key1 in TaxaPair_Reln_Dict:
				#DistMat[j][i] = DistMat[i][j] = TaxaPair_Reln_Dict[key1]._GetXLSumGeneTrees()
				DistMat[j][i] = DistMat[i][j] = TaxaPair_Reln_Dict[key1]._GetNormalizedXLSumGeneTrees()
			elif key2 in TaxaPair_Reln_Dict:
				#DistMat[j][i] = DistMat[i][j] = TaxaPair_Reln_Dict[key2]._GetXLSumGeneTrees()
				DistMat[j][i] = DistMat[i][j] = TaxaPair_Reln_Dict[key2]._GetNormalizedXLSumGeneTrees()
			else:
				DistMat[j][i] = DistMat[i][j] = 0
	# end add - sourya
	#---------------------------------------
	# comment - sourya
	## add - sourya
	#for i in range(no_of_clust - 1):
		#for j in range(i+1, no_of_clust):
			## here both clust_species_list[i] and clust_species_list[j]
			## are one element lists (according to their construction)
			## we have extracted the corresponding element by using [0] operator (extracting first element)
			#curr_taxa_pair_list = []
			#for k1 in range(len(clust_species_list[i])):
	#for k2 in range(len(clust_species_list[j])):	  
		#x1 = clust_species_list[i][k1]
		#x2 = clust_species_list[j][k2]
		#key1 = (x1, x2)
		#key2 = (x2, x1)
		##print 'key1: ', key1, ' key2: ', key2
		#if key1 in TaxaPair_Reln_Dict:
			##DistMat[j][i] = DistMat[i][j] = TaxaPair_Reln_Dict[key1]._GetXLSumGeneTrees()
				#curr_taxa_pair_list.append(TaxaPair_Reln_Dict[key1]._GetNormalizedXLSumGeneTrees())
		#else:
			##DistMat[j][i] = DistMat[i][j] = TaxaPair_Reln_Dict[key2]._GetXLSumGeneTrees()
			#curr_taxa_pair_list.append(TaxaPair_Reln_Dict[key2]._GetNormalizedXLSumGeneTrees())
			## average of this pairwise list is used as the XT approximation
			#DistMat[j][i] = DistMat[i][j] = (sum(curr_taxa_pair_list) * 1.0) / len(curr_taxa_pair_list)
	## end add - sourya
	# end comment - sourya
	#---------------------------------------
					
	# loop to execute the agglomerative clustering
	while(no_of_clust > 2): 
				
		if (DEBUG_LEVEL > 2):
			fp = open(Output_Text_File, 'a')
			fp.write('\n iteration start --- number of clusters: ' + str(no_of_clust))
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
				# we normalize the extra lineage based score
				# by the sum of extra lineages for all other taxa from the taxa indices i and j
				Norm_DistMat[i][j] = ObtainNormalizedVal(DistMat[i][j], sum_list[i], sum_list[j])
				Norm_DistMat[j][i] = Norm_DistMat[i][j]
		
		if (DEBUG_LEVEL > 2):
			fp = open(Output_Text_File, 'a')
			fp.write('\n printing contents of sum_list --- ' + str(sum_list))
			fp.close()
			PrintMatrixContent(no_of_clust, clust_species_list, Norm_DistMat, 'Norm_DistMat', Output_Text_File)
		
		# add - sourya
		# now we have to find the minimum among these elements 
		# present in the matrix Norm_DistMat
		min_val = Norm_DistMat[0][1]
		min_idx_i = 0
		min_idx_j = 1
		for i in range(no_of_clust - 1):
			for j in range(i+1, no_of_clust):
				if (i == j):
					continue
				if (Norm_DistMat[i][j] < min_val):
					min_val = Norm_DistMat[i][j]
					min_idx_i = i
					min_idx_j = j
				elif (Norm_DistMat[i][j] == min_val):
					# here we prioritize the cluster pair having minimum number of species
					if (len(clust_species_list[i]) + len(clust_species_list[j])) < (len(clust_species_list[min_idx_i]) + len(clust_species_list[min_idx_j])):
						min_idx_i = i
						min_idx_j = j
		# end add - sourya

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
			Curr_tree.update_splits(delete_outdegree_one=False)
			
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
			Curr_tree.update_splits(delete_outdegree_one=False)
			
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
			Curr_tree.update_splits(delete_outdegree_one=False)
			
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
			Curr_tree.update_splits(delete_outdegree_one=False)
		
		#---------------------------------------------------------      
		if (DEBUG_LEVEL > 2):
			fp = open(Output_Text_File, 'a')
			fp.write('\n label of newnode: ' + str(Node_Label(newnode)))
			fp.write('\n label of all taxa mrca node (recomputed): ' + str(Node_Label(Curr_tree.mrca(taxon_labels=taxa_list))))      
			fp.write('\n before inserting row col, DistMat dimension: ' + str(DistMat.size))
			fp.close()
		
		# adjust the DistMat by inserting one new row and column corresponding to the new cluster
		# and then deleting the information of earlier two clusters
		# first append one row
		DistMat = numpy.vstack((DistMat, numpy.zeros((1, no_of_clust), dtype=numpy.float)))
		# then append one column
		DistMat = numpy.hstack((DistMat, numpy.zeros((no_of_clust + 1, 1), dtype=numpy.float)))
		
		if (DEBUG_LEVEL > 2):
			fp = open(Output_Text_File, 'a')
			fp.write('\n after inserting row col, DistMat dimension: ' + str(DistMat.size))
			fp.close()
			
		DistMat = numpy.reshape(DistMat, ((no_of_clust + 1), (no_of_clust + 1)), order='C')
		
		# MRCA based adjustment
		# now fill the elements of the new added row and column
		for k in range(no_of_clust):
			# for any index k, the number of extra lineage is the maximum of these three quantities
			DistMat[k][no_of_clust] = max(DistMat[k][min_idx_i], DistMat[k][min_idx_j], DistMat[min_idx_i][min_idx_j])
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
def Refine_Supertree_Binary_Form(Gene_TreeList, Curr_tree, Output_Text_File):
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
			ResolveMultifurcation(Gene_TreeList, Curr_tree, clust_species_list, len(curr_node_children), Output_Text_File)


