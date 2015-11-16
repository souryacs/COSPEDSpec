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
		# get the out edge list of the current node which are not explored yet 
		outnodes = []
		
		# add - sourya
		# if all the descendants of current cluster are leaves then we turn the flag on
		flag_all_leaf_desc = True
		# end add - sourya
		
		for l in Cluster_Info_Dict[root_clust_node_idx]._GetOutEdgeList():
			if (Cluster_Info_Dict[l]._GetExploredStatus() == 0):
				outnodes.append(l)
				# add - sourya
				if (len(Cluster_Info_Dict[l]._GetOutEdgeList()) > 0):
					# this cluster is non leaf
					flag_all_leaf_desc = False
				# end add - sourya
			
		if (len(outnodes) == 0):
			spec_list = Cluster_Info_Dict[root_clust_node_idx]._GetSpeciesList()
			if (len(spec_list) > 1):
				Tree_Str_List = Tree_Str_List + '('
			Tree_Str_List = Tree_Str_List + ','.join("'" + item + "'" for item in spec_list)
			if (len(spec_list) > 1):
				Tree_Str_List = Tree_Str_List + ')'
		else:
			
			# comment - sourya
			#Tree_Str_List = Tree_Str_List + '('
			
			Tree_Str_List = Tree_Str_List + ','.join("'" + item + "'" for item in Cluster_Info_Dict[root_clust_node_idx]._GetSpeciesList())
			Tree_Str_List = Tree_Str_List + ','    
			
			if (flag_all_leaf_desc == False):	# this condition add - sourya
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
			
			if (flag_all_leaf_desc == False):	# this condition add - sourya
				Tree_Str_List = Tree_Str_List + ')'
				
			# comment - sourya
			#Tree_Str_List = Tree_Str_List + ')'
		
	return Tree_Str_List    

#--------------------------------------------------------
# this function defines relationship between a pair of nodes in a tree
# the relationship is either ancestor / descendant, or siblings, or no relationship 
def DefineLeafPairReln(xl_val, node1, node2, reln_type, Curr_tree_taxa_count):	# extra parameter add - sourya
	node1_level = node1.level()
	node2_level = node2.level()
	
	key1 = (node1.taxon.label, node2.taxon.label)
	key2 = (node2.taxon.label, node1.taxon.label)
	
	if key1 in TaxaPair_Reln_Dict:
		TaxaPair_Reln_Dict[key1]._AddSupportingTree()
		TaxaPair_Reln_Dict[key1]._AddXLVal(xl_val)
		TaxaPair_Reln_Dict[key1]._AddEdgeCount(reln_type)
		if (node1_level < node2_level):
			#if (reln_type == RELATION_R4):
				#TaxaPair_Reln_Dict[key1]._IncrLevelDiffInfoCount(0, (node2_level - node1_level))
			# comment - sourya
			#TaxaPair_Reln_Dict[key1]._IncrAllRelnLevelDiffInfoCount(0, (node2_level - node1_level))
			# add - sourya
			TaxaPair_Reln_Dict[key1]._IncrAllRelnLevelDiffInfoCount(0, ((node2_level - node1_level) * 1.0) / Curr_tree_taxa_count)
			# end add - sourya
		elif (node1_level > node2_level):
			#if (reln_type == RELATION_R4):
				#TaxaPair_Reln_Dict[key1]._IncrLevelDiffInfoCount(1, (node1_level - node2_level))
			# comment - sourya
			#TaxaPair_Reln_Dict[key1]._IncrAllRelnLevelDiffInfoCount(1, (node1_level - node2_level))
			# add - sourya
			TaxaPair_Reln_Dict[key1]._IncrAllRelnLevelDiffInfoCount(1, ((node1_level - node2_level) * 1.0) / Curr_tree_taxa_count)
			# end add - sourya
		else:
			#if (reln_type == RELATION_R4):
				#TaxaPair_Reln_Dict[key1]._IncrLevelDiffInfoCount(2, 0)
			TaxaPair_Reln_Dict[key1]._IncrAllRelnLevelDiffInfoCount(2, 0)
	elif key2 in TaxaPair_Reln_Dict:
		TaxaPair_Reln_Dict[key2]._AddSupportingTree()
		TaxaPair_Reln_Dict[key2]._AddXLVal(xl_val)
		TaxaPair_Reln_Dict[key2]._AddEdgeCount(Complementary_Reln(reln_type))
		if (node1_level < node2_level):
			#if (reln_type == RELATION_R4):
				#TaxaPair_Reln_Dict[key2]._IncrLevelDiffInfoCount(1, (node2_level - node1_level))
			# comment - sourya
			#TaxaPair_Reln_Dict[key2]._IncrAllRelnLevelDiffInfoCount(1, (node2_level - node1_level))
			# add - sourya
			TaxaPair_Reln_Dict[key2]._IncrAllRelnLevelDiffInfoCount(1, ((node2_level - node1_level) * 1.0) / Curr_tree_taxa_count)
			# end add - sourya
		elif (node1_level > node2_level):
			#if (reln_type == RELATION_R4):
				#TaxaPair_Reln_Dict[key2]._IncrLevelDiffInfoCount(0, (node1_level - node2_level))
			# comment - sourya
			#TaxaPair_Reln_Dict[key2]._IncrAllRelnLevelDiffInfoCount(0, (node1_level - node2_level))
			# add - sourya
			TaxaPair_Reln_Dict[key2]._IncrAllRelnLevelDiffInfoCount(0, ((node1_level - node2_level) * 1.0) / Curr_tree_taxa_count)
			# end add - sourya
		else:
			#if (reln_type == RELATION_R4):
				#TaxaPair_Reln_Dict[key2]._IncrLevelDiffInfoCount(2, 0)
			TaxaPair_Reln_Dict[key2]._IncrAllRelnLevelDiffInfoCount(2, 0)
	else:
		TaxaPair_Reln_Dict.setdefault(key1, Reln_TaxaPair())
		TaxaPair_Reln_Dict[key1]._AddSupportingTree()
		TaxaPair_Reln_Dict[key1]._AddXLVal(xl_val)
		TaxaPair_Reln_Dict[key1]._AddEdgeCount(reln_type)
		if (node1_level < node2_level):
			#if (reln_type == RELATION_R4):
				#TaxaPair_Reln_Dict[key1]._IncrLevelDiffInfoCount(0, (node2_level - node1_level))
			# comment - sourya
			#TaxaPair_Reln_Dict[key1]._IncrAllRelnLevelDiffInfoCount(0, (node2_level - node1_level))
			# add - sourya
			TaxaPair_Reln_Dict[key1]._IncrAllRelnLevelDiffInfoCount(0, ((node2_level - node1_level) * 1.0) / Curr_tree_taxa_count)
			# end add - sourya
		elif (node1_level > node2_level):
			#if (reln_type == RELATION_R4):
				#TaxaPair_Reln_Dict[key1]._IncrLevelDiffInfoCount(1, (node1_level - node2_level))
			# comment - sourya
			#TaxaPair_Reln_Dict[key1]._IncrAllRelnLevelDiffInfoCount(1, (node1_level - node2_level))
			# add - sourya
			TaxaPair_Reln_Dict[key1]._IncrAllRelnLevelDiffInfoCount(1, ((node1_level - node2_level) * 1.0) / Curr_tree_taxa_count)
			# end add - sourya
		else:
			#if (reln_type == RELATION_R4):
				#TaxaPair_Reln_Dict[key1]._IncrLevelDiffInfoCount(2, 0)
			TaxaPair_Reln_Dict[key1]._IncrAllRelnLevelDiffInfoCount(2, 0)
			
	return

#--------------------------------------------------------
# this function derives coupket relations belonging to one tree
# that is provided as an input argument to this function
def DeriveCoupletRelations(Curr_tree):
  
	Curr_tree_taxa_count = len(Curr_tree.infer_taxa().labels())

	# traverse the internal nodes of the tree in postorder fashion
	for curr_node in Curr_tree.postorder_internal_node_iter():        
		# this is the level value associated with this node
		#curr_node_level = curr_node.level()
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
					DefineLeafPairReln(xl_val, curr_node_child_leaf_nodes[i], curr_node_child_leaf_nodes[j], RELATION_R3, Curr_tree_taxa_count)
		
		# one leaf node (direct descendant) and another leaf node (under one internal node)
		# will be related by ancestor / descendant relations
		if (len(curr_node_child_leaf_nodes) > 0) and (len(curr_node_child_internal_nodes) > 0):
			for p in curr_node_child_leaf_nodes:
				for q in curr_node_child_internal_nodes:
					for r in q.leaf_nodes():
						DefineLeafPairReln(xl_val, p, r, RELATION_R1, Curr_tree_taxa_count)
		
		# finally a pair of leaf nodes which are descendant of internal nodes will be related by RELATION_R4 relation
		if (len(curr_node_child_internal_nodes) > 1):
			for i in range(len(curr_node_child_internal_nodes) - 1):
				for j in range(i+1, len(curr_node_child_internal_nodes)):
					for p in curr_node_child_internal_nodes[i].leaf_nodes():
						for q in curr_node_child_internal_nodes[j].leaf_nodes():
							DefineLeafPairReln(xl_val, p, q, RELATION_R4, Curr_tree_taxa_count)

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

#-----------------------------------------------------
# this is the taxa list generated from current internal node
def GetTaxaUnderInternalNode(curr_node):
	taxa_list_from_curr_internal_node = []
	for n in curr_node.leaf_nodes():
		taxa_list_from_curr_internal_node.append(n.taxon.label)
	return taxa_list_from_curr_internal_node

#----------------------------------------
def Complementary_Reln(inp_reln):
  if (inp_reln == RELATION_R3) or (inp_reln == RELATION_R4):
    return inp_reln
  elif (inp_reln == RELATION_R1):
    return RELATION_R2
  else:
    return RELATION_R1
