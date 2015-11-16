import Header
from Header import *

"""
this document contains functions to manage the taxa clusters
regarding their connectivity, cost, and update operations
"""

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
				score_val = score_val + TaxaPair_Reln_Dict[key1]._GetConnPrVal(RELATION_R1)
				# add - sourya
				# now we employ frequency based scoring mechanism, to satisfy maximum agreement property
				#score_val = score_val + TaxaPair_Reln_Dict[key1]._GetEdgeWeight(RELATION_R1)
				# end add - sourya
			elif key2 in TaxaPair_Reln_Dict:
				# comment - sourya
				# previously the priority metric based scoring was employed to determine the scoring 
				# among two taxa clusters	
				score_val = score_val + TaxaPair_Reln_Dict[key2]._GetConnPrVal(RELATION_R2)
				# add - sourya
				# now we employ frequency based scoring mechanism, to satisfy maximum agreement property
				#score_val = score_val + TaxaPair_Reln_Dict[key2]._GetEdgeWeight(RELATION_R2)
				# end add - sourya
			else:
				if (DEBUG_LEVEL > 2):
					fp.write('\n score compute -- key pair ' + str(t1) + ',' + str(t2) + ' does not exist ')

	if (DEBUG_LEVEL > 2):
		fp.write('\n pairwise score of this cluster pair is : ' + str(score_val))
		fp.close()
		
	return score_val


#-----------------------------------------------------    
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

#-----------------------------------------------------        
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

#-----------------------------------------------------  
''' 
this function performs transitive reduction of a graph (transitive closure) and subsequently modifies the cluster of nodes
in terms of the edge connectivity, to make it free of redunant edges 
'''
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
    
#-----------------------------------------------------
""" this function creates one new cluster with the given index value
also, it inserts one specified taxa in that cluster """
def Create_Cluster_Taxa_Label(target_clust_idx, target_taxa_label):
	# create the cluster
	Cluster_Info_Dict.setdefault(target_clust_idx, Cluster_node(target_taxa_label))
	# include the cluster idx in the global list CURRENT_CLUST_IDX_LIST
	CURRENT_CLUST_IDX_LIST.append(target_clust_idx)
	# mention the cluster index in the taxa information
	Taxa_Info_Dict[target_taxa_label]._Set_Clust_Idx_taxa_Part(target_clust_idx)

#-----------------------------------------------------
""" this function appends one specified taxon on a given cluster """
def Append_Cluster_Taxa_Label(target_clust_idx, target_taxa_label):
	if target_taxa_label not in Cluster_Info_Dict[target_clust_idx]._GetSpeciesList():
		Cluster_Info_Dict[target_clust_idx]._Append_taxa(target_taxa_label)
		# mention the cluster index in the taxa information
		Taxa_Info_Dict[target_taxa_label]._Set_Clust_Idx_taxa_Part(target_clust_idx)  
	