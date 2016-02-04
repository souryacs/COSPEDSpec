"""
this document contains functions to manage the taxa clusters
regarding their connectivity, cost, and update operations
"""

import Header
from Header import *
import UtilFunc
from UtilFunc import *

#-----------------------------------------------------
"""
this function generates an out edge from the first cluster to the second cluster
"""
def ConnectClustPairOutEdge(Reachability_Graph_Mat, clust1, clust2):
	Cluster_Info_Dict[clust1]._AddOutEdge(clust2)
	Cluster_Info_Dict[clust2]._AddInEdge(clust1)
	reach_clust1_idx = CURRENT_CLUST_IDX_LIST.index(clust1)
	reach_clust2_idx = CURRENT_CLUST_IDX_LIST.index(clust2)
	Reachability_Graph_Mat[reach_clust1_idx][reach_clust2_idx] = 1
	return

#-----------------------------------------------------
"""
this function finds whether there is a hidden R1 / R2 relation between this pair of clusters
which was originally related via R4 (no edge)
the return value is an integer. It can have following values.
@ 0: no edge will be established
@ 1: there will be a directed out edge from the clust1_key to the clust2_key
@ 2: there will be a directed out edge from the clust2_key to the clust1_key
@ 3: there will be a directed out edge from the parent(s) of clust1_key to the clust2_key
@ 4: there will be a directed out edge from the parent(s) of clust2_key to the clust1_key
"""
def FindClusterReln(clust1_key, clust2_key):
	"""
	now there is no such strict consensus R4 relation between this cluster pair
	so we can proceed to check if a edge connection can be possible
	"""
	clust1_spec_list = Cluster_Info_Dict[clust1_key]._GetSpeciesList()
	clust2_spec_list = Cluster_Info_Dict[clust2_key]._GetSpeciesList()

	"""
	there are two conditions for a cluster to have directed out edge to another cluster
	-- For a couplet (x,y) where x in Clust1, y in Clust2
	1) (Either fr1(x,y) + pr1(x,y) - pr2(x,y)) > (fr4(x,y) - pr1(x,y))
	or 2) Level ratio of r1(x,y) is greater than a certain threshold, and fr1(x,y) is significant (r1 relation 
	belongs to the set of allowable relations from x to y)
	"""
	
	for t1 in clust1_spec_list:
		for t2 in clust2_spec_list:
			key1 = (t1, t2)
			key2 = (t2, t1)
			if key1 in TaxaPair_Reln_Dict:
				r1_freq = TaxaPair_Reln_Dict[key1]._GetEdgeWeight(RELATION_R1)
				r2_freq = TaxaPair_Reln_Dict[key1]._GetEdgeWeight(RELATION_R2)
				r4_freq = TaxaPair_Reln_Dict[key1]._GetEdgeWeight(RELATION_R4)
				pseudo_r1_freq = TaxaPair_Reln_Dict[key1]._GetFreqPseudoR1(0)
				pseudo_r2_freq = TaxaPair_Reln_Dict[key1]._GetFreqPseudoR1(1)
				r1_level_val_ratio = TaxaPair_Reln_Dict[key1]._GetLevelValRatio(0)
				r2_level_val_ratio = TaxaPair_Reln_Dict[key1]._GetLevelValRatio(1)
				allowed_reln_list = TaxaPair_Reln_Dict[key1]._GetAllowedRelnList()

				if ((r1_freq + 2 * (pseudo_r1_freq - pseudo_r2_freq)) >= r4_freq) \
					and ((round(r1_level_val_ratio, 2) >= R1R2Reln_MAJ_THRS_low) and (RELATION_R1 in allowed_reln_list)):
					#if (len(clust1_spec_list) == 1) and (len(clust2_spec_list) == 1):
						#if (round(((r1_freq * 1.0) / TaxaPair_Reln_Dict[key1]._GetConsensusFreq()), 2) >= CONSENSUS_FREQ_RATIO_THR) \
							#or (round(r1_level_val_ratio, 2) >= R1R2Reln_MAJ_THRS_high):
							#return 1
					if (len(clust1_spec_list) > 1):
						return 1
					else:
						return 3

				if ((r2_freq + 2 * (pseudo_r2_freq - pseudo_r1_freq)) >= r4_freq) \
					and ((round(r2_level_val_ratio, 2) >= R1R2Reln_MAJ_THRS_low) and (RELATION_R2 in allowed_reln_list)):
					#if (len(clust1_spec_list) == 1) and (len(clust2_spec_list) == 1):
						#if (round(((r2_freq * 1.0) / TaxaPair_Reln_Dict[key1]._GetConsensusFreq()), 2) >= CONSENSUS_FREQ_RATIO_THR) \
							#or (round(r2_level_val_ratio, 2) >= R1R2Reln_MAJ_THRS_high):
							#return 2
					if (len(clust2_spec_list) > 1):
						return 2
					else:
						return 4

			elif key2 in TaxaPair_Reln_Dict:
				r1_freq = TaxaPair_Reln_Dict[key2]._GetEdgeWeight(RELATION_R1)
				r2_freq = TaxaPair_Reln_Dict[key2]._GetEdgeWeight(RELATION_R2)
				r4_freq = TaxaPair_Reln_Dict[key2]._GetEdgeWeight(RELATION_R4)
				pseudo_r1_freq = TaxaPair_Reln_Dict[key2]._GetFreqPseudoR1(0)
				pseudo_r2_freq = TaxaPair_Reln_Dict[key2]._GetFreqPseudoR1(1)
				r1_level_val_ratio = TaxaPair_Reln_Dict[key2]._GetLevelValRatio(0)
				r2_level_val_ratio = TaxaPair_Reln_Dict[key2]._GetLevelValRatio(1)
				allowed_reln_list = TaxaPair_Reln_Dict[key2]._GetAllowedRelnList()
				
				if ((r1_freq + 2 * (pseudo_r1_freq - pseudo_r2_freq)) >= r4_freq) \
					and ((round(r1_level_val_ratio, 2) >= R1R2Reln_MAJ_THRS_low) and (RELATION_R1 in allowed_reln_list)):
					#if (len(clust1_spec_list) == 1) and (len(clust2_spec_list) == 1):
						#if (round(((r1_freq * 1.0) / TaxaPair_Reln_Dict[key2]._GetConsensusFreq()), 2) >= CONSENSUS_FREQ_RATIO_THR) \
							#or (round(r1_level_val_ratio, 2) >= R1R2Reln_MAJ_THRS_high):
							#return 2
					if (len(clust2_spec_list) > 1):
						return 2
					else:
						return 4

				if ((r2_freq + 2 * (pseudo_r2_freq - pseudo_r1_freq)) >= r4_freq) \
					and ((round(r2_level_val_ratio, 2) >= R1R2Reln_MAJ_THRS_low) and (RELATION_R2 in allowed_reln_list)):
					#if (len(clust1_spec_list) == 1) and (len(clust2_spec_list) == 1):
						#if (round(((r2_freq * 1.0) / TaxaPair_Reln_Dict[key2]._GetConsensusFreq()), 2) >= CONSENSUS_FREQ_RATIO_THR) \
							#or (round(r2_level_val_ratio, 2) >= R1R2Reln_MAJ_THRS_high):
							#return 1
					if (len(clust1_spec_list) > 1):
						return 1
					else:
						return 3

	return 0

#-----------------------------------------------------
# this function computes the score (ancestor relation) from clust1 to clust2
def ComputeScore(clust1, clust2, Output_Text_File, MPP_SOLVE_METRIC, DIST_MAT_TYPE):
  
	if (DEBUG_LEVEL > 2):
		fp = open(Output_Text_File, 'a')
		fp.write('\n ===>>> Score --- cluster 1 - taxa set: ' + str(Cluster_Info_Dict[clust1]._GetSpeciesList())\
			+ ' cluster 2 taxa set ' + str(Cluster_Info_Dict[clust2]._GetSpeciesList()))
		
	score_val = 0
	couplet_count = 0
	
	for t1 in Cluster_Info_Dict[clust1]._GetSpeciesList():
		for t2 in Cluster_Info_Dict[clust2]._GetSpeciesList():
			key1 = (t1, t2)
			key2 = (t2, t1)
			if key1 in TaxaPair_Reln_Dict:
				couplet_count = couplet_count + 1
				# modified - sourya - using normalization by the number of support trees
				if (MPP_SOLVE_METRIC == 1):
					score_val = score_val + (TaxaPair_Reln_Dict[key1]._GetConnPrVal(RELATION_R1) * 1.0) / TaxaPair_Reln_Dict[key1]._GetNoSupportTrees()
				else:
					score_val = score_val + TaxaPair_Reln_Dict[key1]._GetNormalizedXLSumGeneTrees(DIST_MAT_TYPE)
			elif key2 in TaxaPair_Reln_Dict:
				couplet_count = couplet_count + 1
				# modified - sourya - using normalization by the number of support trees
				if (MPP_SOLVE_METRIC == 1):
					score_val = score_val + (TaxaPair_Reln_Dict[key2]._GetConnPrVal(RELATION_R2) * 1.0) / TaxaPair_Reln_Dict[key2]._GetNoSupportTrees()
				else:
					score_val = score_val + TaxaPair_Reln_Dict[key2]._GetNormalizedXLSumGeneTrees(DIST_MAT_TYPE)
			else:
				if (DEBUG_LEVEL > 2):
					fp.write('\n score compute -- key pair ' + str(t1) + ',' + str(t2) + ' does not exist ')

	if (DEBUG_LEVEL > 2):
		fp.write('\n pairwise score of this cluster pair is : ' + str((score_val * 1.0) / couplet_count))
		fp.close()
		
	return (score_val * 1.0) / couplet_count


#-----------------------------------------------------    
""" this function solves multiple parenting problem (C2)
by uniquely selecting one particular parent
the selection is carried out using a scoring mechanism """
def SolveMultipleParentC2Problem(Output_Text_File, MPP_SOLVE_METRIC, DIST_MAT_TYPE):
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
				scoring_dict[cz] = ComputeScore(cz, cx, Output_Text_File, MPP_SOLVE_METRIC, DIST_MAT_TYPE)
			
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

			if (MPP_SOLVE_METRIC == 1):
				# for the priority measure 
				# after sorting the scoring list, now remove all except the last element from the 
				# cx cluster in edge lists
				for i in range(len(Scoring_List) - 1):
					target_delete_clust_idx = Scoring_List[i][0]
					Cluster_Info_Dict[cx]._RemoveInEdge(target_delete_clust_idx)
					Cluster_Info_Dict[target_delete_clust_idx]._RemoveOutEdge(cx)
			else:
				# for the XL based measure
				# after sorting the scoring list, now remove all except the first element from the 
				# cx cluster in edge lists
				for i in range(1, len(Scoring_List)):
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
			# A->B case
			if (Reachability_Graph_Mat[i][j] == 1):
				for k in range(no_of_clusters):
					# A->C and B->C case
					if (Reachability_Graph_Mat[j][k] == 1) and (Reachability_Graph_Mat[i][k] == 1):
						# comment - sourya - check
						#Reachability_Graph_Mat[i][k] = 0
						
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
