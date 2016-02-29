#!/usr/bin/env python

"""
this file contains code for detection of conflicts generated due to possible inclusion of one new relation
in the supertree 
the conflict can be generated due to mismatch of the relationships
or due to the mismatch of the branch length settings 
"""

import Header
from Header import *

#-------------------------------------------------------
"""
this function checks whether the given couplet has a connection already
in terms of the Reachability_Graph_Mat entries
which are filled by the transitive connection
"""
def CheckExistingConn(clust1, clust2, Reachability_Graph_Mat, target_reln, Output_Text_File):
	if (target_reln == RELATION_R2):
		src_taxa_clust_idx = clust2
		dest_taxa_clust_idx = clust1
	else:
		src_taxa_clust_idx = clust1
		dest_taxa_clust_idx = clust2

	# check whether the taxa pairs are already resolved
		
	# if the taxa pair is already part of the same cluster
	# that is, they are already related via equivalence partition 
	# then return
	if (src_taxa_clust_idx == dest_taxa_clust_idx):
		if (DEBUG_LEVEL >= 2):
			fp = open(Output_Text_File, 'a')
			fp.write('\n target_reln: ' + str(target_reln) + ' already in same taxa cluster')
			fp.close()
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
				fp.write('\n src cluster: ' + str(src_taxa_clust_idx) + '  dest cluster: ' + str(dest_taxa_clust_idx) + ' target_reln: ' + str(target_reln) + ' already related via relation r1')
			elif (Reachability_Graph_Mat[reach_mat_src_taxa_clust_idx][reach_mat_dest_taxa_clust_idx] == 2):
				fp.write('\n src cluster: ' + str(src_taxa_clust_idx) + '  dest cluster: ' + str(dest_taxa_clust_idx) + ' target_reln: ' + str(target_reln) + ' already related via relation r4')
			elif (Reachability_Graph_Mat[reach_mat_dest_taxa_clust_idx][reach_mat_src_taxa_clust_idx] == 1):
				fp.write('\n src cluster: ' + str(src_taxa_clust_idx) + '  dest cluster: ' + str(dest_taxa_clust_idx) + ' target_reln: ' + str(target_reln) + ' already related via relation r2')
			fp.close()
		return 1

	return 0

#-------------------------------------------------------
"""
this function checks whether a given relation between a couplet 
induces transitive conflict (with respect to other relations between other couplets
which are already included in the set of included relations S_r)
"""
def CheckTransitiveConflict(clust1, clust2, Reachability_Graph_Mat, target_reln, Output_Text_File):
	if (target_reln == RELATION_R2):
		src_taxa_clust_idx = clust2
		dest_taxa_clust_idx = clust1
	else:
		src_taxa_clust_idx = clust1
		dest_taxa_clust_idx = clust2
	
	# we find the indices of the reachability matrix corresponding to individual cluster indices 
	reach_mat_src_taxa_clust_idx = CURRENT_CLUST_IDX_LIST.index(src_taxa_clust_idx)
	reach_mat_dest_taxa_clust_idx = CURRENT_CLUST_IDX_LIST.index(dest_taxa_clust_idx)
	
	# check for the transitive conflict
	if (target_reln == RELATION_R2) or (target_reln == RELATION_R1):
		"""
		if A->B is to be established
		and there exists D->A 
		then if B->D then return a conflict
		"""
		for x in Cluster_Info_Dict[src_taxa_clust_idx]._GetInEdgeList():
			if (Reachability_Graph_Mat[reach_mat_dest_taxa_clust_idx][CURRENT_CLUST_IDX_LIST.index(x)] == 1):
				if (DEBUG_LEVEL >= 2):
					fp = open(Output_Text_File, 'a')
					fp.write('\n Conflict  - condition 1')
					fp.close()
				return 1
		
		"""
		if A->B is to be established
		and there exists B->E 
		then if E->A then return a conflict 
		"""
		for x in Cluster_Info_Dict[dest_taxa_clust_idx]._GetOutEdgeList():
			if (Reachability_Graph_Mat[CURRENT_CLUST_IDX_LIST.index(x)][reach_mat_src_taxa_clust_idx] == 1):
				if (DEBUG_LEVEL >= 2):
					fp = open(Output_Text_File, 'a')
					fp.write('\n Conflict  - condition 2')
					fp.close()
				return 1
	
		"""
		if A->B is to be established
		and there exists D->A and B->E 
		then if E->D or E=D then return a conflict  
		"""
		for x in Cluster_Info_Dict[src_taxa_clust_idx]._GetInEdgeList():
			for y in Cluster_Info_Dict[dest_taxa_clust_idx]._GetOutEdgeList():
				if (x == y) or (Reachability_Graph_Mat[CURRENT_CLUST_IDX_LIST.index(y)][CURRENT_CLUST_IDX_LIST.index(x)] == 1):
					if (DEBUG_LEVEL >= 2):
						fp = open(Output_Text_File, 'a')
						fp.write('\n Conflict  - condition 3')
						fp.close()
					return 1
	
		"""
		if A->B is to be established
		and there exists D><A
		if there is B->D or D=B then return a conflict  
		"""
		for x in Cluster_Info_Dict[src_taxa_clust_idx]._GetNoEdgeList():
			if (dest_taxa_clust_idx == x):	# or (Reachability_Graph_Mat[reach_mat_dest_taxa_clust_idx][CURRENT_CLUST_IDX_LIST.index(x)] == 1):
				if (DEBUG_LEVEL >= 2):
					fp = open(Output_Text_File, 'a')
					fp.write('\n Conflict  - condition 4')
					fp.close()
				return 1
	
		#"""
		#if A->B is to be established
		#and there exists D><A
		#then for all B->E
		#so if E->D or D=E then return a conflict  
		#"""
		#for x in Cluster_Info_Dict[src_taxa_clust_idx]._GetNoEdgeList():
			#for y in Cluster_Info_Dict[dest_taxa_clust_idx]._GetOutEdgeList():
				#if (x == y) or (Reachability_Graph_Mat[CURRENT_CLUST_IDX_LIST.index(y)][CURRENT_CLUST_IDX_LIST.index(x)] == 1):
					#if (DEBUG_LEVEL >= 2):
						#fp = open(Output_Text_File, 'a')
						#fp.write('\n Conflict  - condition 5')
						#fp.close()
					#return 1
	
	elif (target_reln == RELATION_R3):
		"""
		# if A=B is to be established
		# and there exists B->D 
		# then if D->A or D=A then return a conflict    
		"""
		for x in Cluster_Info_Dict[dest_taxa_clust_idx]._GetOutEdgeList():
			if (x == src_taxa_clust_idx) or (Reachability_Graph_Mat[CURRENT_CLUST_IDX_LIST.index(x)][reach_mat_src_taxa_clust_idx] == 1):
				if (DEBUG_LEVEL >= 2):
					fp = open(Output_Text_File, 'a')
					fp.write('\n Conflict  - condition 6')
					fp.close()
				return 1
			
		"""
		if A=B is to be established
		and there exists A->D 
		then if D->B or D=B then return a conflict    	
		"""
		for x in Cluster_Info_Dict[src_taxa_clust_idx]._GetOutEdgeList():
			if (x == dest_taxa_clust_idx) or (Reachability_Graph_Mat[CURRENT_CLUST_IDX_LIST.index(x)][reach_mat_dest_taxa_clust_idx] == 1):
				if (DEBUG_LEVEL >= 2):
					fp = open(Output_Text_File, 'a')
					fp.write('\n Conflict  - condition 7')
					fp.close()
				return 1
		
		"""
		# if A=B is to be established
		# and there exists D->B 
		# then if A->D or A=D then return a conflict    	
		"""
		for x in Cluster_Info_Dict[dest_taxa_clust_idx]._GetInEdgeList():
			if (x == src_taxa_clust_idx) or (Reachability_Graph_Mat[reach_mat_src_taxa_clust_idx][CURRENT_CLUST_IDX_LIST.index(x)] == 1):
				if (DEBUG_LEVEL >= 2):
					fp = open(Output_Text_File, 'a')
					fp.write('\n Conflict  - condition 8')
					fp.close()
				return 1

		"""
		if A=B is to be established
		and there exists D->A 
		then if B->D or B=D then return a conflict    	
		"""
		for x in Cluster_Info_Dict[src_taxa_clust_idx]._GetInEdgeList():
			if (x == dest_taxa_clust_idx) or (Reachability_Graph_Mat[reach_mat_dest_taxa_clust_idx][CURRENT_CLUST_IDX_LIST.index(x)] == 1):
				if (DEBUG_LEVEL >= 2):
					fp = open(Output_Text_File, 'a')
					fp.write('\n Conflict  - condition 9')
					fp.close()
				return 1
		
		"""
		if A=B is to be established
		and there exists B><D 
		then if A=D then return a conflict    	
		"""
		for x in Cluster_Info_Dict[dest_taxa_clust_idx]._GetNoEdgeList():
			if (x == src_taxa_clust_idx):
				if (DEBUG_LEVEL >= 2):
					fp = open(Output_Text_File, 'a')
					fp.write('\n Conflict  - condition 10')
					fp.close()
				return 1
		"""
		if A=B is to be established
		and there exists A><D 
		then if B=D then return a conflict    	
		"""
		for x in Cluster_Info_Dict[src_taxa_clust_idx]._GetNoEdgeList():
			if (x == dest_taxa_clust_idx):
				if (DEBUG_LEVEL >= 2):
					fp = open(Output_Text_File, 'a')
					fp.write('\n Conflict  - condition 11')
					fp.close()
				return 1
	
	return 0

