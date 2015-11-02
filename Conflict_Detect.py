#!/usr/bin/env python

"""
this file contains code for detection of conflicts generated due to possible inclusion of one new relation
in the supertree 
the conflict can be generated due to mismatch of the relationships
or due to the mismatch of the branch length settings 
"""

import Header
from Header import *
import Cost_Update
from Cost_Update import *
import ReachGraph_Update
from ReachGraph_Update import *
import UtilFunc
from UtilFunc import *

#-------------------------------------------------------
""" 
this function checks whether the current relation (given as an input) is feasible to the input supertree configuration
only the current relations in the supertree is used for conflict detection 
"""
def Possible_Conflict_Curr_Reln(src_taxa_label, dest_taxa_label, Reachability_Graph_Mat, target_reln, Output_Text_File):
	if (target_reln == RELATION_R2):
		src_taxa_clust_idx = Taxa_Info_Dict[dest_taxa_label]._Get_Taxa_Part_Clust_Idx()
		dest_taxa_clust_idx = Taxa_Info_Dict[src_taxa_label]._Get_Taxa_Part_Clust_Idx()
	else:
		src_taxa_clust_idx = Taxa_Info_Dict[src_taxa_label]._Get_Taxa_Part_Clust_Idx()
		dest_taxa_clust_idx = Taxa_Info_Dict[dest_taxa_label]._Get_Taxa_Part_Clust_Idx()

	#--------------------------------------------------------------------  
	# check whether the taxa pairs are already resolved
		
	# if the taxa pair is already part of the same cluster
	# that is, they are already related via equivalence partition 
	# then return
	if (src_taxa_clust_idx == dest_taxa_clust_idx):
		return 1
		
	# we find the indices of the reachability matrix corresponding to individual cluster indices 
	reach_mat_src_taxa_clust_idx = CURRENT_CLUST_IDX_LIST.index(src_taxa_clust_idx)
	reach_mat_dest_taxa_clust_idx = CURRENT_CLUST_IDX_LIST.index(dest_taxa_clust_idx)

	# case A - if the clusters are already related (depicted in the reachability matrix) then return 1
	if (Reachability_Graph_Mat[reach_mat_dest_taxa_clust_idx][reach_mat_src_taxa_clust_idx] > 0) or \
			(Reachability_Graph_Mat[reach_mat_src_taxa_clust_idx][reach_mat_dest_taxa_clust_idx] > 0):
		if (DEBUG_LEVEL >= 2):
			fp = open(Output_Text_File, 'a')
			if (Reachability_Graph_Mat[reach_mat_src_taxa_clust_idx][reach_mat_dest_taxa_clust_idx] == 1):
				fp.write('\n target_reln: ' + str(target_reln) + ' already related via relation r1')
			elif (Reachability_Graph_Mat[reach_mat_src_taxa_clust_idx][reach_mat_dest_taxa_clust_idx] == 2):
				fp.write('\n target_reln: ' + str(target_reln) + ' already related via relation r4')
			elif (Reachability_Graph_Mat[reach_mat_dest_taxa_clust_idx][reach_mat_src_taxa_clust_idx] == 1):
				fp.write('\n target_reln: ' + str(target_reln) + ' already related via relation r2')
			fp.close()
		return 1

	#--------------------------------------------------------------------
	# check for the transitive conflict
	if (target_reln == RELATION_R2) or (target_reln == RELATION_R1):
		# if A->B is to be established
		# and there exists D->A 
		# then if B->D or B><D then return a conflict
		for x in Cluster_Info_Dict[src_taxa_clust_idx]._GetInEdgeList():
			if (Reachability_Graph_Mat[reach_mat_dest_taxa_clust_idx][CURRENT_CLUST_IDX_LIST.index(x)] > 0):
				if (DEBUG_LEVEL >= 2):
					fp = open(Output_Text_File, 'a')
					fp.write('\n **** CASE 1 *** ')
					fp.close()
				return 1
		
		# if A->B is to be established
		# and there exists B->E 
		# then if E->A or E><A then return a conflict  
		for x in Cluster_Info_Dict[dest_taxa_clust_idx]._GetOutEdgeList():
			if (Reachability_Graph_Mat[CURRENT_CLUST_IDX_LIST.index(x)][reach_mat_src_taxa_clust_idx] > 0):
				if (DEBUG_LEVEL >= 2):
					fp = open(Output_Text_File, 'a')
					fp.write('\n **** CASE 2 *** ')
					fp.close()
				return 1

		# if A->B is to be established
		# and there exists D->A and B->E 
		# then if E->D or E=D or E><D then return a conflict  
		for x in Cluster_Info_Dict[src_taxa_clust_idx]._GetInEdgeList():
			for y in Cluster_Info_Dict[dest_taxa_clust_idx]._GetOutEdgeList():
				if (x == y) or (Reachability_Graph_Mat[CURRENT_CLUST_IDX_LIST.index(y)][CURRENT_CLUST_IDX_LIST.index(x)] > 0):
					if (DEBUG_LEVEL >= 2):
						fp = open(Output_Text_File, 'a')
						fp.write('\n **** CASE 3 *** ')
						fp.close()
					return 1

		# if A->B is to be established
		# and there exists D><A
		# then it would be D><B afterwards 
		# so if there is D->B, B->D or D=B then return a conflict  
		for x in Cluster_Info_Dict[src_taxa_clust_idx]._GetNoEdgeList():
			if (dest_taxa_clust_idx == x) \
				or (Reachability_Graph_Mat[CURRENT_CLUST_IDX_LIST.index(x)][reach_mat_dest_taxa_clust_idx] == 1) \
				or (Reachability_Graph_Mat[reach_mat_dest_taxa_clust_idx][CURRENT_CLUST_IDX_LIST.index(x)] == 1):
				if (DEBUG_LEVEL >= 2):
					fp = open(Output_Text_File, 'a')
					fp.write('\n **** CASE 4 *** ')
					fp.close()
				return 1
		
		# if A->B is to be established
		# and there exists D><A
		# then for all B->E
		# it would be D><E afterwards 
		# so if D->E or E->D or D=E then return a conflict  
		for x in Cluster_Info_Dict[src_taxa_clust_idx]._GetNoEdgeList():
			for y in Cluster_Info_Dict[dest_taxa_clust_idx]._GetOutEdgeList():
				if (x == y) \
					or (Reachability_Graph_Mat[CURRENT_CLUST_IDX_LIST.index(y)][CURRENT_CLUST_IDX_LIST.index(x)] == 1) \
					or (Reachability_Graph_Mat[CURRENT_CLUST_IDX_LIST.index(x)][CURRENT_CLUST_IDX_LIST.index(y)] == 1):
					if (DEBUG_LEVEL >= 2):
						fp = open(Output_Text_File, 'a')
						fp.write('\n **** CASE 5 *** ')
						fp.close()
					return 1

	elif (target_reln == RELATION_R3):
		# if A=B is to be established
		# and there exists B->D 
		# then if D->A or D><A or D=A then return a conflict    
		for x in Cluster_Info_Dict[dest_taxa_clust_idx]._GetOutEdgeList():
			if (x == src_taxa_clust_idx) or (Reachability_Graph_Mat[CURRENT_CLUST_IDX_LIST.index(x)][reach_mat_src_taxa_clust_idx] > 0):
				return 1
		# if A=B is to be established
		# and there exists D->B 
		# then if A->D or A><D or A=D then return a conflict    	
		for x in Cluster_Info_Dict[dest_taxa_clust_idx]._GetInEdgeList():
			if (x == src_taxa_clust_idx) or (Reachability_Graph_Mat[reach_mat_src_taxa_clust_idx][CURRENT_CLUST_IDX_LIST.index(x)] > 0):
				return 1
		
		# if A=B is to be established
		# and there exists B><D 
		# then if A->D or D->A or A=D then return a conflict    	
		for x in Cluster_Info_Dict[dest_taxa_clust_idx]._GetNoEdgeList():
			if (x == src_taxa_clust_idx) or (Reachability_Graph_Mat[reach_mat_src_taxa_clust_idx][CURRENT_CLUST_IDX_LIST.index(x)] == 1) \
				or (Reachability_Graph_Mat[CURRENT_CLUST_IDX_LIST.index(x)][reach_mat_src_taxa_clust_idx] == 1):
				return 1

		# if A=B is to be established
		# and there exists A->D 
		# then if D->B or D><B or D=B then return a conflict    	
		for x in Cluster_Info_Dict[src_taxa_clust_idx]._GetOutEdgeList():
			if (x == dest_taxa_clust_idx) or (Reachability_Graph_Mat[CURRENT_CLUST_IDX_LIST.index(x)][reach_mat_dest_taxa_clust_idx] > 0):
				return 1
		# if A=B is to be established
		# and there exists D->A 
		# then if B->D or B><D or B=D then return a conflict    	
		for x in Cluster_Info_Dict[src_taxa_clust_idx]._GetInEdgeList():
			if (x == dest_taxa_clust_idx) or (Reachability_Graph_Mat[reach_mat_dest_taxa_clust_idx][CURRENT_CLUST_IDX_LIST.index(x)] > 0):
				return 1
		
		# if A=B is to be established
		# and there exists A><D 
		# then if D->B or B->D or B=D then return a conflict    	
		for x in Cluster_Info_Dict[src_taxa_clust_idx]._GetNoEdgeList():
			if (x == dest_taxa_clust_idx) or (Reachability_Graph_Mat[reach_mat_dest_taxa_clust_idx][CURRENT_CLUST_IDX_LIST.index(x)] == 1) \
				or (Reachability_Graph_Mat[CURRENT_CLUST_IDX_LIST.index(x)][reach_mat_dest_taxa_clust_idx] == 1):
				return 1
		
	else:	#if (target_reln == RELATION_R4):
		# construct the out neighborhood of src_cluster
		# it will contain the cluster itself and all the other clusters connected via out edges from this cluster
		src_clust_out_neighb = []
		src_clust_out_neighb.append(src_taxa_clust_idx)
		src_clust_out_neighb.extend(Cluster_Info_Dict[src_taxa_clust_idx]._GetOutEdgeList())
		# construct the out neighborhood of dest cluster
		# it will contain the cluster itself and all the other clusters connected via out edges from this cluster    
		dest_clust_out_neighb = []
		dest_clust_out_neighb.append(dest_taxa_clust_idx)
		dest_clust_out_neighb.extend(Cluster_Info_Dict[dest_taxa_clust_idx]._GetOutEdgeList())
		
		# if mutually any pair of taxa belonging to respective out edge clusters
		# are themselves related via any relationships other than NO EDGE then this connection is not possible
		for x in src_clust_out_neighb:
			for y in dest_clust_out_neighb:
				# comment - sourya
				#if (x == y) or (Reachability_Graph_Mat[CURRENT_CLUST_IDX_LIST.index(y)][CURRENT_CLUST_IDX_LIST.index(x)] == 1) \
					#or (Reachability_Graph_Mat[CURRENT_CLUST_IDX_LIST.index(x)][CURRENT_CLUST_IDX_LIST.index(y)] == 1):
					#return 1
				# add - sourya
				if (x == y):
					if (DEBUG_LEVEL >= 2):
						fp = open(Output_Text_File, 'a')
						fp.write('\n **** Same cluster --- with species list : ' + str(Cluster_Info_Dict[x]._GetSpeciesList()))
						fp.close()
					return 1
				if (Reachability_Graph_Mat[CURRENT_CLUST_IDX_LIST.index(y)][CURRENT_CLUST_IDX_LIST.index(x)] == 1):
					if (DEBUG_LEVEL >= 2):
						fp = open(Output_Text_File, 'a')
						fp.write('\n **** Out edge from the cluster ( ' + str(Cluster_Info_Dict[y]._GetSpeciesList()) \
							+ ' ) to the cluster ( ' + str(Cluster_Info_Dict[x]._GetSpeciesList()) + ' )')
						fp.close()
					return 1
				if (Reachability_Graph_Mat[CURRENT_CLUST_IDX_LIST.index(x)][CURRENT_CLUST_IDX_LIST.index(y)] == 1):
					if (DEBUG_LEVEL >= 2):
						fp = open(Output_Text_File, 'a')
						fp.write('\n **** Out edge from the cluster ( ' + str(Cluster_Info_Dict[x]._GetSpeciesList()) \
							+ ' ) to the cluster ( ' + str(Cluster_Info_Dict[y]._GetSpeciesList()) + ' )')
						fp.close()
					return 1
				# end add - sourya
	
	return 0    
    