#!/usr/bin/env python

import Header
from Header import *    
  
##-----------------------------------------------------
## original function - sourya
##---------------------------------
## this function prints the tree in Newick format
## sourya - this is the old function with cluster based species list printing in the newick format
#def PrintNewick_Original(root_clust_node_idx):
	#if 0:
		#print 'in function printnewick:   root_clust_node_idx: ', root_clust_node_idx
		#print 'taxa set: ', Cluster_Info_Dict[root_clust_node_idx]._GetSpeciesList()  
		#print 'out clust list: ', Cluster_Info_Dict[root_clust_node_idx]._GetOutEdgeList()

	#Tree_Str_List = ''
	## process the node provided it has not been explored yet
	#if (Cluster_Info_Dict[root_clust_node_idx]._GetExploredStatus() == 0):  
		## set the explored status of the current node to true
		#Cluster_Info_Dict[root_clust_node_idx]._SetExploredStatus()
		## get the out edge list of the current node which are not explored yet 
		#outnodes = []
		#for l in Cluster_Info_Dict[root_clust_node_idx]._GetOutEdgeList():
			#if (Cluster_Info_Dict[l]._GetExploredStatus() == 0):
				#outnodes.append(l)
		#if (len(outnodes) == 0):
			#spec_list = Cluster_Info_Dict[root_clust_node_idx]._GetSpeciesList()
			#if (len(spec_list) > 1):
				#Tree_Str_List = Tree_Str_List + '('
			#Tree_Str_List = Tree_Str_List + ','.join("'" + item + "'" for item in spec_list)
			#if (len(spec_list) > 1):
				#Tree_Str_List = Tree_Str_List + ')'
		#else:
			#Tree_Str_List = Tree_Str_List + '('
			#Tree_Str_List = Tree_Str_List + ','.join("'" + item + "'" for item in Cluster_Info_Dict[root_clust_node_idx]._GetSpeciesList())
			#Tree_Str_List = Tree_Str_List + ','    
			#Tree_Str_List = Tree_Str_List + '('
			#for i in range(len(outnodes)):
				#if (Cluster_Info_Dict[outnodes[i]]._GetExploredStatus() == 0):  
					#Tree_Str_List = Tree_Str_List + PrintNewick(outnodes[i])
					#if (i < (len(outnodes) - 1)):
						## we check whether any subsequent node belonging to the outnodes list
						## is left for traverse
						#j = i + 1
						#while (j < len(outnodes)):
							#if (Cluster_Info_Dict[outnodes[j]]._GetExploredStatus() == 0):  
								#break
							#j = j + 1	      
						## in this case, we append one comma
						#if (j < len(outnodes)):
							#Tree_Str_List = Tree_Str_List + ','      
			
			#Tree_Str_List = Tree_Str_List + ')'
			#Tree_Str_List = Tree_Str_List + ')'
		
	#return Tree_Str_List    
    
## end original function - sourya
#-------------------------------------------

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
		# this is the list of taxa of this cluster
		spec_list = Cluster_Info_Dict[root_clust_node_idx]._GetSpeciesList()
		# get the out edge list of the current cluster which are not explored yet 
		outnodes = []
		for l in Cluster_Info_Dict[root_clust_node_idx]._GetOutEdgeList():
			if (Cluster_Info_Dict[l]._GetExploredStatus() == 0):
				outnodes.append(l)
		
		""" 
		at first, print the contents of this taxa cluster
		if the cluster has more than one taxon, then use ( and ) to enclose the taxa list
		"""
		if (len(outnodes) > 0):	# and (len(spec_list) == 1):
			Tree_Str_List = Tree_Str_List + '('
			
		if (len(spec_list) > 1):
			Tree_Str_List = Tree_Str_List + '('
		Tree_Str_List = Tree_Str_List + ','.join("'" + item + "'" for item in spec_list)
		if (len(spec_list) > 1):
			Tree_Str_List = Tree_Str_List + ')'
		"""
		here we check if the cluster has one or more out edges
		then recursively traverse all the out edge clusters
		"""
		if (len(outnodes) > 0):
			# first add one comma
			Tree_Str_List = Tree_Str_List + ','
			
			# then add one opening bracket, within which, all the out edge cluster contents will reside
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
			
			# at last, append one closing bracket, signifying the end of out edge cluster contents
			Tree_Str_List = Tree_Str_List + ')'

		if (len(outnodes) > 0):	# and (len(spec_list) == 1):
			Tree_Str_List = Tree_Str_List + ')'
		
		
		
		## add - sourya
		#single_tax_single_leaf_outnode = False
		## end add - sourya
		
		#if (len(outnodes) == 0):
			#if (len(spec_list) > 1):
				#Tree_Str_List = Tree_Str_List + '('
			#Tree_Str_List = Tree_Str_List + ','.join("'" + item + "'" for item in spec_list)
			#if (len(spec_list) > 1):
				#Tree_Str_List = Tree_Str_List + ')'
		#else:
			
			## add - sourya
			#if (len(spec_list) == 1):
				#if (len(outnodes) == 1):
					#l = outnodes[0]
					#if (len(Cluster_Info_Dict[l]._GetSpeciesList()) == 1) and (len(Cluster_Info_Dict[l]._GetOutEdgeList()) == 0):
						#single_tax_single_leaf_outnode = True
			## end add - sourya
			
			### add - sourya
			##if (len(spec_list) > 1):
				##single_tax_single_leaf_outnode = False
			##elif (len(outnodes) > 1):
				##single_tax_single_leaf_outnode = False
			##else:
				##l = outnodes[0]
				##if (len(Cluster_Info_Dict[l]._GetOutEdgeList()) > 0):
					##single_tax_single_leaf_outnode = False
			### end add - sourya
			
			## comment - sourya
			#if (single_tax_single_leaf_outnode == False):	#1:	#(len(spec_list) > 1):	# condition add - sourya
				#Tree_Str_List = Tree_Str_List + '('
			## end comment - sourya
			
			#Tree_Str_List = Tree_Str_List + ','.join("'" + item + "'" for item in spec_list)
			#Tree_Str_List = Tree_Str_List + ','    
			
			#if (single_tax_single_leaf_outnode == False):	# 1:	#(flag_all_leaf_desc == False):	# this condition add - sourya
				#Tree_Str_List = Tree_Str_List + '('
			
			#for i in range(len(outnodes)):
				#if (Cluster_Info_Dict[outnodes[i]]._GetExploredStatus() == 0):  
					#Tree_Str_List = Tree_Str_List + PrintNewick(outnodes[i])
					#if (i < (len(outnodes) - 1)):
						## we check whether any subsequent node belonging to the outnodes list
						## is left for traverse
						#j = i + 1
						#while (j < len(outnodes)):
							#if (Cluster_Info_Dict[outnodes[j]]._GetExploredStatus() == 0):  
								#break
							#j = j + 1
						## in this case, we append one comma
						#if (j < len(outnodes)):
							#Tree_Str_List = Tree_Str_List + ','      
			
			#if (single_tax_single_leaf_outnode == False):	#1:	#(flag_all_leaf_desc == False):	# this condition add - sourya
				#Tree_Str_List = Tree_Str_List + ')'
				
			## comment - sourya
			#if (single_tax_single_leaf_outnode == False):	#1:	#(len(spec_list) > 1):	# condition add - sourya
				#Tree_Str_List = Tree_Str_List + ')'
			## end comment - sourya
		
		
	return Tree_Str_List    

#--------------------------------------------------------
# this function defines relationship between a pair of nodes in a tree
# the relationship is either ancestor / descendant, or siblings, or no relationship 
def DefineLeafPairReln(xl_val, ratio_val, lca_level, node1, node2, reln_type, Curr_tree_taxa_count):	# extra parameter add - sourya
	node1_level = node1.level()
	node2_level = node2.level()
	
	key1 = (node1.taxon.label, node2.taxon.label)
	key2 = (node2.taxon.label, node1.taxon.label)
	
	if key1 in TaxaPair_Reln_Dict:
		TaxaPair_Reln_Dict[key1]._AddSupportingTree()
		TaxaPair_Reln_Dict[key1]._AddXLVal(xl_val)
		TaxaPair_Reln_Dict[key1]._AddEdgeCount(reln_type, ratio_val)	#sourya
		# modified - sourya
		if (node1_level < node2_level):
			TaxaPair_Reln_Dict[key1]._IncrAllRelnLevelDiffInfoCount(0, ((node2_level - node1_level) * 1.0) / Curr_tree_taxa_count)
		elif (node1_level > node2_level):
			TaxaPair_Reln_Dict[key1]._IncrAllRelnLevelDiffInfoCount(1, ((node1_level - node2_level) * 1.0) / Curr_tree_taxa_count)
		else:	#if (node1_level == node2_level):
			TaxaPair_Reln_Dict[key1]._IncrAllRelnLevelDiffInfoCount(2, 0)
		#-----------------------
		# add - sourya
		if (reln_type == RELATION_R4):
			if ((node1_level - lca_level) == 2) and ((node2_level - lca_level) > 2):
				TaxaPair_Reln_Dict[key1]._AddFreqPseudoR1(0, ratio_val)	#sourya
			elif ((node1_level - lca_level) > 2) and ((node2_level - lca_level) == 2):
				TaxaPair_Reln_Dict[key1]._AddFreqPseudoR1(1, ratio_val)	#sourya
		#-----------------------
	elif key2 in TaxaPair_Reln_Dict:
		TaxaPair_Reln_Dict[key2]._AddSupportingTree()
		TaxaPair_Reln_Dict[key2]._AddXLVal(xl_val)
		TaxaPair_Reln_Dict[key2]._AddEdgeCount(Complementary_Reln(reln_type), ratio_val)	#sourya
		# modified - sourya
		if (node1_level < node2_level):
			TaxaPair_Reln_Dict[key2]._IncrAllRelnLevelDiffInfoCount(1, ((node2_level - node1_level) * 1.0) / Curr_tree_taxa_count)
		elif (node1_level > node2_level):
			TaxaPair_Reln_Dict[key2]._IncrAllRelnLevelDiffInfoCount(0, ((node1_level - node2_level) * 1.0) / Curr_tree_taxa_count)
		else:	#if (node1_level == node2_level):
			TaxaPair_Reln_Dict[key2]._IncrAllRelnLevelDiffInfoCount(2, 0)
		#-----------------------
		# add - sourya
		if (reln_type == RELATION_R4):
			if ((node1_level - lca_level) == 2) and ((node2_level - lca_level) > 2):
				TaxaPair_Reln_Dict[key2]._AddFreqPseudoR1(1, ratio_val)	#sourya
			elif ((node1_level - lca_level) > 2) and ((node2_level - lca_level) == 2):
				TaxaPair_Reln_Dict[key2]._AddFreqPseudoR1(0, ratio_val)	#sourya
		#-----------------------
	else:
		TaxaPair_Reln_Dict.setdefault(key1, Reln_TaxaPair())
		TaxaPair_Reln_Dict[key1]._AddSupportingTree()
		TaxaPair_Reln_Dict[key1]._AddXLVal(xl_val)
		TaxaPair_Reln_Dict[key1]._AddEdgeCount(reln_type, ratio_val)	#sourya
		# modified - sourya
		if (node1_level < node2_level):
			TaxaPair_Reln_Dict[key1]._IncrAllRelnLevelDiffInfoCount(0, ((node2_level - node1_level) * 1.0) / Curr_tree_taxa_count)
		elif (node1_level > node2_level):
			TaxaPair_Reln_Dict[key1]._IncrAllRelnLevelDiffInfoCount(1, ((node1_level - node2_level) * 1.0) / Curr_tree_taxa_count)
		else:	#if (node1_level == node2_level):
			TaxaPair_Reln_Dict[key1]._IncrAllRelnLevelDiffInfoCount(2, 0)
		#-----------------------
		# add - sourya
		if (reln_type == RELATION_R4):
			if ((node1_level - lca_level) == 2) and ((node2_level - lca_level) > 2):
				TaxaPair_Reln_Dict[key1]._AddFreqPseudoR1(0, ratio_val)	#sourya
			elif ((node1_level - lca_level) > 2) and ((node2_level - lca_level) == 2):
				TaxaPair_Reln_Dict[key1]._AddFreqPseudoR1(1, ratio_val)	#sourya
		#-----------------------
			
	return

#--------------------------------------------------------
# this function derives coupket relations belonging to one tree
# that is provided as an input argument to this function
def DeriveCoupletRelations(Curr_tree, Total_Taxa_Count):
  
	Curr_tree_taxa_count = len(Curr_tree.infer_taxa().labels())
	
	# modified  - sourya - 
	# initially it was introduced, but again reverted back due to debugging
	ratio_val = 1
	#ratio_val = (Curr_tree_taxa_count * 1.0) / Total_Taxa_Count

	# traverse the internal nodes of the tree in postorder fashion
	for curr_node in Curr_tree.postorder_internal_node_iter():        
		# this is the level value associated with this node
		curr_node_level = curr_node.level()
		# compute the XL value associated with this node
		#xl_val = (len(curr_node.leaf_nodes()) - 2)
		xl_val = ((len(curr_node.leaf_nodes()) - 2) * 1.0 ) / Curr_tree_taxa_count
		
		# list the leaf and internal children of the current node
		curr_node_child_leaf_nodes = []
		curr_node_child_internal_nodes = []
		for x in curr_node.child_nodes():
			if (x.is_leaf() == True):
				curr_node_child_leaf_nodes.append(x)
			else:
				curr_node_child_internal_nodes.append(x)
		
		# pair of leaf nodes will be related by sibling relations
		if (len(curr_node_child_leaf_nodes) > 1):
			for i in range(len(curr_node_child_leaf_nodes) - 1):
				for j in range(i+1, len(curr_node_child_leaf_nodes)):
					DefineLeafPairReln(xl_val, ratio_val, curr_node_level, curr_node_child_leaf_nodes[i], curr_node_child_leaf_nodes[j], \
						RELATION_R3, Curr_tree_taxa_count)
		
		# one leaf node (direct descendant) and another leaf node (under one internal node)
		# will be related by ancestor / descendant relations
		if (len(curr_node_child_leaf_nodes) > 0) and (len(curr_node_child_internal_nodes) > 0):
			for p in curr_node_child_leaf_nodes:
				for q in curr_node_child_internal_nodes:
					for r in q.leaf_nodes():
						DefineLeafPairReln(xl_val, ratio_val, curr_node_level, p, r, RELATION_R1, Curr_tree_taxa_count)
		
		# finally a pair of leaf nodes which are descendant of internal nodes will be related by RELATION_R4 relation
		if (len(curr_node_child_internal_nodes) > 1):
			for i in range(len(curr_node_child_internal_nodes) - 1):
				for j in range(i+1, len(curr_node_child_internal_nodes)):
					for p in curr_node_child_internal_nodes[i].leaf_nodes():
						for q in curr_node_child_internal_nodes[j].leaf_nodes():
							DefineLeafPairReln(xl_val, ratio_val, curr_node_level, p, q, RELATION_R4, Curr_tree_taxa_count)

##-----------------------------------------------------
# this function reads the input tree list file
# parameters: ROOTED_TREE - whether the treelist to be read as rooted format
# PRESERVE_UNDERSCORE: whether underscores of the taxa name will be preserved or not
# INPUT_FILE_FORMAT: data is read from the file according to NEWICK or NEXUS format
# INPUT_FILENAME: file containing the input treelist

def Read_Input_Treelist(ROOTED_TREE, PRESERVE_UNDERSCORE, INPUT_FILE_FORMAT, INPUT_FILENAME):
	Inp_TreeList = dendropy.TreeList.get_from_path(INPUT_FILENAME, schema=INPUT_FILE_FORMAT, \
							preserve_underscores=PRESERVE_UNDERSCORE, \
							default_as_rooted=ROOTED_TREE)

	return Inp_TreeList

##-----------------------------------------------------
# this function reads an input tree from a specified file
# parameters: ROOTED_TREE - whether the treelist to be read as rooted format
# PRESERVE_UNDERSCORE: whether underscores of the taxa name will be preserved or not
# INPUT_FILE_FORMAT: data is read from the file according to NEWICK or NEXUS format
# INPUT_FILENAME: file containing the input tree

def Read_Input_Tree(ROOTED_TREE, PRESERVE_UNDERSCORE, INPUT_FILE_FORMAT, INPUT_FILENAME):
	Inp_Tree = dendropy.Tree.get_from_path(INPUT_FILENAME, schema=INPUT_FILE_FORMAT, \
							preserve_underscores=PRESERVE_UNDERSCORE, \
							default_as_rooted=ROOTED_TREE)

	return Inp_Tree

##-----------------------------------------------------
# this function writes the treelist in the specified file
# parameters: Inp_TreeList - the treelist to be written
# outfile: target output text file name
# FILE_FORMAT: data is written according to NEWICK (1) or NEXUS (2) format
# Suppress_Root: condition whether the treelist will have rooting information
# Unquote_Underscore: condition whether the treelist will have unquoted underscores

def Write_Output_Treelist(Inp_TreeList, outfile, FILE_FORMAT, Suppress_Root=False, Unquote_Underscore=False):
	Inp_TreeList.write_to_path(outfile, FILE_FORMAT, suppress_rooting=Suppress_Root, unquoted_underscores=Unquote_Underscore)
  
##-----------------------------------------------------
# this function writes a single tree in a specified file
# parameters: Inp_Tree - the tree to be written
# outfile: target output text file name
# FILE_FORMAT: data is written according to NEWICK (1) or NEXUS (2) format
# Suppress_Root: condition whether the tree will have rooting information
# Unquote_Underscore: condition whether the treelist will have unquoted underscores

def Write_Output_Tree(Inp_Tree, outfile, FILE_FORMAT, Suppress_Root=False, Unquote_Underscore=False):
	Inp_Tree.write(open(outfile, 'w'), FILE_FORMAT, suppress_rooting=Suppress_Root, unquoted_underscores=Unquote_Underscore)

#--------------------------------------------------
# this function returns the label of an internal or a leaf node 
# in terms of newick representation
def Node_Label(inp_node):
	return str(inp_node.as_newick_string(suppress_edge_lengths=True))

##-----------------------------------------------------
## this is the taxa list generated from current internal node
#def GetTaxaUnderInternalNode(curr_node):
	#taxa_list_from_curr_internal_node = []
	#for n in curr_node.leaf_nodes():
		#taxa_list_from_curr_internal_node.append(n.taxon.label)
	#return taxa_list_from_curr_internal_node

#----------------------------------------
def Complementary_Reln(inp_reln):
  if (inp_reln == RELATION_R3) or (inp_reln == RELATION_R4):
    return inp_reln
  elif (inp_reln == RELATION_R1):
    return RELATION_R2
  else:
    return RELATION_R1

#------------------------------------------------
"""
this function computes average XL information between a pair of taxa clusters
@param type_of_output: if 0, computes the average of XL measures
													1, returns the minimum of XL measures
													2, returns the maximum of XL measures
"""
def FindAvgXL(taxa_clust1, taxa_clust2, DIST_MAT_TYPE, single_elem, type_of_output=0):
	
	#print '***** taxa_clust1: ', taxa_clust1, ' taxa_clust2: ', taxa_clust2
	
	if (single_elem == False):
		curr_taxa_pair_list = []
	
	for x1 in taxa_clust1:
		for x2 in taxa_clust2:  
			key1 = (x1, x2)
			key2 = (x2, x1)
			#print 'key1: ', key1, ' key2: ', key2
			if key1 in TaxaPair_Reln_Dict:
				val = TaxaPair_Reln_Dict[key1]._GetNormalizedXLSumGeneTrees(DIST_MAT_TYPE)
				if (single_elem == False):
					curr_taxa_pair_list.append(val)
				else:
					return val
			elif key2 in TaxaPair_Reln_Dict:
				val = TaxaPair_Reln_Dict[key2]._GetNormalizedXLSumGeneTrees(DIST_MAT_TYPE)
				if (single_elem == False):
					curr_taxa_pair_list.append(val)
				else:
					return val
	
	#print 'curr_taxa_pair_list : ', curr_taxa_pair_list
	#print 'Avg of curr_taxa_pair_list: ', (sum(curr_taxa_pair_list) * 1.0) / len(curr_taxa_pair_list)
	#print 'Stdev of curr_taxa_pair_list: ', (numpy.std(numpy.array(curr_taxa_pair_list)))
	
	# average of this pairwise list is used as the XL approximation
	if (single_elem == False):
		if (len(curr_taxa_pair_list) > 0):
			if (type_of_output == 0):
				return (sum(curr_taxa_pair_list) * 1.0) / len(curr_taxa_pair_list)
			elif (type_of_output == 1):
				return min(curr_taxa_pair_list)
			else:
				return max(curr_taxa_pair_list)
		else:
			return 0
	
	return 0
