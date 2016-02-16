#!/usr/bin/env python

import Header
from Header import *
import Cost_Update
from Cost_Update import *
import ReachGraph_Update
from ReachGraph_Update import *
import Conflict_Detect
from Conflict_Detect import *

#--------------------------------------------------------
"""
this function processes the max priority queue
containing support score values for different couplets
here, only positive support scores are processed
also, establishment of only R1 / R2 relations are considered
"""
def Proc_Queue_Pos_Score_R1R2(Reachability_Graph_Mat, Output_Text_File):
	Inp_Queue = Cost_List_Taxa_Pair_Multi_Reln
	while (0 < len(Inp_Queue)):
		""" 
		extract the 1st element of "Inp_Queue" 
		since it is sorted to have max support score value at the beginning 
		"""
		outlist = Heap_Extract_Max(Inp_Queue)
		
		src_taxa_label = outlist[0]
		dest_taxa_label = outlist[1]
		reln_type = outlist[2]
		conn_score = TaxaPair_Reln_Dict[(src_taxa_label, dest_taxa_label)]._GetEdgeCost_ConnReln(reln_type)
		
		if (DEBUG_LEVEL >= 2):
			fp = open(Output_Text_File, 'a')
			fp.write('\n ===>> CONFLICTING QUEUE -- ')      
			fp.write(' current extracted max element: ' + str(src_taxa_label) + ' and ' + str(dest_taxa_label) + \
								' reln type: ' + str(reln_type) + ' conn score: ' + str(conn_score))
			fp.close()

		"""
		R3 relations are already processed
		"""
		if (reln_type == RELATION_R3):
			continue

		"""
		Case A - relation is either R1 or R2
		and the relation is consensus
		in such a case, immediately apply the relation
		provided the relation is not conflicting
		"""
		if ((reln_type == RELATION_R1) or (reln_type == RELATION_R2)):
			"""
			check whether it is a consensus relation
			"""
			if (TaxaPair_Reln_Dict[(src_taxa_label, dest_taxa_label)]._CheckTargetRelnConsensus(reln_type) == True):
				""" 
				if the couplet is already connected, as shown by the entries in Reachability_Graph_Mat
				then there is no need for any connection
				"""
				existing_conn = CheckExistingConn(src_taxa_label, dest_taxa_label, Reachability_Graph_Mat, reln_type, Output_Text_File)
				if (existing_conn == 0):
					"""
					there is no apparent existing relationship between the couplet
					"""
					"""
					first, check whether the relation induces conflict
					otherwise, we apply the relations
					"""
					if (CheckTransitiveConflict(src_taxa_label, dest_taxa_label, Reachability_Graph_Mat, reln_type, Output_Text_File) == 0):
						if (DEBUG_LEVEL >= 2):
							fp = open(Output_Text_File, 'a')    
							fp.write('\n ==>>>>>>>>> NEW CONN --- nodes to be connected: ' \
								+ str(src_taxa_label) + ' and ' + str(dest_taxa_label) + \
										' nodes indices: ' + str(COMPLETE_INPUT_TAXA_LIST.index(src_taxa_label)) \
											+ ' and ' + str(COMPLETE_INPUT_TAXA_LIST.index(dest_taxa_label)) + \
										' relation type: ' + str(reln_type) + ' conn score: ' + str(conn_score))
							fp.close()
						
						# also update the reachability graph information
						Reachability_Graph_Mat = AdjustReachGraph(Reachability_Graph_Mat, src_taxa_label, dest_taxa_label, reln_type)
				
		"""
		case B - relation is R4 and it is also consensus score
		"""
		if (reln_type == RELATION_R4):
			"""
			check whether it is a consensus relation
			"""
			if (TaxaPair_Reln_Dict[(src_taxa_label, dest_taxa_label)]._CheckTargetRelnConsensus(reln_type) == True):
				""" 
				if the couplet is already connected, as shown by the entries in Reachability_Graph_Mat
				then there is no need for any connection
				"""
				existing_conn = CheckExistingConn(src_taxa_label, dest_taxa_label, Reachability_Graph_Mat, reln_type, Output_Text_File)
				if (existing_conn == 0):
					"""
					there is no apparent existing relationship between the couplet
					"""
					"""
					here we check about possible R1 or R2 relation between the cluster pairs
					otherwise we do not process the couplet for now
					"""
					conn_done, target_reln = CheckHiddenR1R2Reln(Reachability_Graph_Mat, src_taxa_label, dest_taxa_label, Output_Text_File)
					if (conn_done == 1):
						if (DEBUG_LEVEL >= 2):
							fp = open(Output_Text_File, 'a')    
							fp.write('\n ==>>>>>>>>> NEW CONN --- nodes to be connected: ' \
								+ str(src_taxa_label) + ' and ' + str(dest_taxa_label) + \
								' nodes indices: ' + str(COMPLETE_INPUT_TAXA_LIST.index(src_taxa_label)) \
								+ ' and ' + str(COMPLETE_INPUT_TAXA_LIST.index(dest_taxa_label)) + \
								' relation type: ' + str(target_reln) + ' conn score: ' \
								+ str(TaxaPair_Reln_Dict[(src_taxa_label, dest_taxa_label)]._GetEdgeCost_ConnReln(target_reln)))
							fp.close()

					# also update the reachability graph information
					Reachability_Graph_Mat = AdjustReachGraph(Reachability_Graph_Mat, src_taxa_label, dest_taxa_label, target_reln)

	return Reachability_Graph_Mat

#-------------------------------------------------------
"""
this function checks whether R1 relation from clust1_taxa_list to clust2_taxa_list 
can be established or not, when the consensus relation is R4
"""
def CheckR1Reln(clust1_spec_list, clust2_spec_list):
	for t1 in clust1_spec_list:
		for t2 in clust2_spec_list:
			key1 = (t1, t2)
			key2 = (t2, t1)
			if key1 in TaxaPair_Reln_Dict:
				r1_freq = TaxaPair_Reln_Dict[key1]._GetEdgeWeight(RELATION_R1)
				r4_freq = TaxaPair_Reln_Dict[key1]._GetEdgeWeight(RELATION_R4)
				pseudo_r1_freq = TaxaPair_Reln_Dict[key1]._GetFreqPseudoR1(0)
				pseudo_r2_freq = TaxaPair_Reln_Dict[key1]._GetFreqPseudoR1(1)
				r1_level_val_ratio = TaxaPair_Reln_Dict[key1]._GetLevelValRatio(0)
				allowed_reln_list = TaxaPair_Reln_Dict[key1]._GetAllowedRelnList()
				r1_score = TaxaPair_Reln_Dict[key1]._GetEdgeCost_ConnReln(RELATION_R1)
				r4_score = TaxaPair_Reln_Dict[key1]._GetEdgeCost_ConnReln(RELATION_R4)
				r4_priority = TaxaPair_Reln_Dict[key1]._GetConnPrVal(RELATION_R4)
				
				#if 0:	#(len(clust1_spec_list) == 1) and (len(clust2_spec_list) == 1):
					#if ((r1_score < 0) and (r4_score < 0)) or (r4_priority <= (PRIORITY_PERCENT * r4_freq)):
						## condition 1
						#if ((r1_freq + 2 * (pseudo_r1_freq - pseudo_r2_freq)) >= r4_freq) \
							#and (round(r1_level_val_ratio, 2) >= R1R2Reln_MAJ_THRS_low) and (RELATION_R1 in allowed_reln_list):
							#return 1
						## condition 2
						#if ((r1_freq + pseudo_r1_freq - pseudo_r2_freq) >= r4_freq) and (RELATION_R1 in allowed_reln_list):
							#return 1
						## condition 3
						#if (round(r1_level_val_ratio, 2) >= R1R2Reln_MAJ_THRS_high) and (RELATION_R1 in allowed_reln_list):
							#return 1
						
				#else:
				"""
				at least one of the clusters have cardinality > 1
				"""
				if (((r1_freq + 2 * (pseudo_r1_freq - pseudo_r2_freq)) >= r4_freq) \
					or (TaxaPair_Reln_Dict[key1]._CheckTargetRelnLevelConsensus(RELATION_R1, 1))) \
					and ((round(r1_level_val_ratio, 2) >= R1R2Reln_MAJ_THRS_low) and (RELATION_R1 in allowed_reln_list)):
					return 1
	
			if key2 in TaxaPair_Reln_Dict:
				r2_freq = TaxaPair_Reln_Dict[key2]._GetEdgeWeight(RELATION_R2)
				r4_freq = TaxaPair_Reln_Dict[key2]._GetEdgeWeight(RELATION_R4)
				pseudo_r1_freq = TaxaPair_Reln_Dict[key2]._GetFreqPseudoR1(0)
				pseudo_r2_freq = TaxaPair_Reln_Dict[key2]._GetFreqPseudoR1(1)
				r2_level_val_ratio = TaxaPair_Reln_Dict[key2]._GetLevelValRatio(1)
				allowed_reln_list = TaxaPair_Reln_Dict[key2]._GetAllowedRelnList()
				r2_score = TaxaPair_Reln_Dict[key2]._GetEdgeCost_ConnReln(RELATION_R2)
				r4_score = TaxaPair_Reln_Dict[key2]._GetEdgeCost_ConnReln(RELATION_R4)
				r4_priority = TaxaPair_Reln_Dict[key2]._GetConnPrVal(RELATION_R4)
	
				#if 0:	#(len(clust1_spec_list) == 1) and (len(clust2_spec_list) == 1):
					#if ((r2_score < 0) and (r4_score < 0)) or (r4_priority <= (PRIORITY_PERCENT * r4_freq)):
						## condition 1
						#if ((r2_freq + 2 * (pseudo_r2_freq - pseudo_r1_freq)) >= r4_freq) \
							#and (round(r2_level_val_ratio, 2) >= R1R2Reln_MAJ_THRS_low) and (RELATION_R2 in allowed_reln_list):
							#return 1
						## condition 2
						#if ((r2_freq + pseudo_r2_freq - pseudo_r1_freq) >= r4_freq) and (RELATION_R2 in allowed_reln_list):
							#return 1
						## condition 3
						#if (round(r2_level_val_ratio, 2) >= R1R2Reln_MAJ_THRS_high) and (RELATION_R2 in allowed_reln_list):
							#return 1
						
				#else:
				"""
				at least one of the clusters have cardinality > 1
				"""
				if (((r2_freq + 2 * (pseudo_r2_freq - pseudo_r1_freq)) >= r4_freq) \
					or (TaxaPair_Reln_Dict[key2]._CheckTargetRelnLevelConsensus(RELATION_R2, 1))) \
					and ((round(r2_level_val_ratio, 2) >= R1R2Reln_MAJ_THRS_low) and (RELATION_R2 in allowed_reln_list)):
					return 1
	
	return 0

#-------------------------------------------------------
"""
this function checks for cluster pairs
whether they can be further analyzed for R1 / R2 relation

return: 1 if R1 reln from clust1 to clust2 is possible
				-1 if R2 reln from clust1 to clust2 is possible
				0 if no such relation is possible
"""
def CheckCandidateR1R2Reln(clust1_spec_list, clust2_spec_list):
	for t1 in clust1_spec_list:
		for t2 in clust2_spec_list:
			key1 = (t1, t2)
			key2 = (t2, t1)
			if key1 in TaxaPair_Reln_Dict:
				r1_level_val_ratio = TaxaPair_Reln_Dict[key1]._GetLevelValRatio(0)
				r2_level_val_ratio = TaxaPair_Reln_Dict[key1]._GetLevelValRatio(1)
				allowed_reln_list = TaxaPair_Reln_Dict[key1]._GetAllowedRelnList()
				if (RELATION_R1 in allowed_reln_list):
					if (round(r1_level_val_ratio, 2) >= R1R2Reln_MAJ_THRS_low):
						return 1

				if (RELATION_R2 in allowed_reln_list):
					if (round(r2_level_val_ratio, 2) >= R1R2Reln_MAJ_THRS_low):
						return -1

			if key2 in TaxaPair_Reln_Dict:
				r2_level_val_ratio = TaxaPair_Reln_Dict[key2]._GetLevelValRatio(1)
				r1_level_val_ratio = TaxaPair_Reln_Dict[key2]._GetLevelValRatio(0)
				allowed_reln_list = TaxaPair_Reln_Dict[key2]._GetAllowedRelnList()
				if (RELATION_R2 in allowed_reln_list):
					if (round(r2_level_val_ratio, 2) >= R1R2Reln_MAJ_THRS_low):
						return 1

				if (RELATION_R1 in allowed_reln_list):
					if (round(r1_level_val_ratio, 2) >= R1R2Reln_MAJ_THRS_low):
						return -1
				
	return 0

#-------------------------------------------------------
"""
this function checks whether a cluster pair can be related by R1 or R2 relation
even if R4 relation is predominant among them
provided no conflict is induced
Here the source taxa cluster has cardinality > 1
"""
def CheckHiddenR1R2Reln(Reachability_Graph_Mat, src_taxa_label, dest_taxa_label, Output_Text_File):
	
	src_taxa_clust_idx = Taxa_Info_Dict[src_taxa_label]._Get_Taxa_Part_Clust_Idx()
	dest_taxa_clust_idx = Taxa_Info_Dict[dest_taxa_label]._Get_Taxa_Part_Clust_Idx()
	"""
	find the taxa list of respective clusters
	"""
	src_cluster_taxa_list = Cluster_Info_Dict[src_taxa_clust_idx]._GetSpeciesList()
	dest_cluster_taxa_list = Cluster_Info_Dict[dest_taxa_clust_idx]._GetSpeciesList()

	"""
	case A - src cluster size > 1
	dest cluster size = 1
	R1 relation from source to destination cluster is sought
	"""
	if (len(src_cluster_taxa_list) > 1) and (len(dest_cluster_taxa_list) == 1):
		if (CheckR1Reln(src_cluster_taxa_list, dest_cluster_taxa_list) == 1):
			if (CheckTransitiveConflict(src_taxa_label, dest_taxa_label, Reachability_Graph_Mat, RELATION_R1, Output_Text_File) == 0):
				return 1, RELATION_R1

	"""
	case B - src cluster size = 1
	dest cluster size > 1
	R2 relation from source to destination cluster is sought
	"""
	if (len(src_cluster_taxa_list) == 1) and (len(dest_cluster_taxa_list) > 1):
		if (CheckR1Reln(dest_cluster_taxa_list, src_cluster_taxa_list) == 1):
			if (CheckTransitiveConflict(src_taxa_label, dest_taxa_label, Reachability_Graph_Mat, RELATION_R2, Output_Text_File) == 0):
				return 1, RELATION_R2

	"""
	case C - src cluster size > 1
	dest cluster size > 1
	R1 / R2 relation from source to destination cluster is sought
	"""
	if ((len(src_cluster_taxa_list) > 1) and (len(dest_cluster_taxa_list) > 1)):
		if (CheckR1Reln(src_cluster_taxa_list, dest_cluster_taxa_list) == 1):
			if (CheckTransitiveConflict(src_taxa_label, dest_taxa_label, Reachability_Graph_Mat, RELATION_R1, Output_Text_File) == 0):
				return 1, RELATION_R1

		if (CheckR1Reln(dest_cluster_taxa_list, src_cluster_taxa_list) == 1):
			if (CheckTransitiveConflict(src_taxa_label, dest_taxa_label, Reachability_Graph_Mat, RELATION_R2, Output_Text_File) == 0):
				return 1, RELATION_R2


	"""
	case D - src cluster size = 1 and dest cluster size = 1
	check R1 / R2 relation between this cluster pair
	"""
	if (len(src_cluster_taxa_list) == 1) and (len(dest_cluster_taxa_list) == 1):
		res = CheckCandidateR1R2Reln(src_cluster_taxa_list, dest_cluster_taxa_list)
		if (res == 1):
			Cluster_Info_Dict[src_taxa_clust_idx]._AddPossibleR1(dest_taxa_clust_idx)
			if src_taxa_clust_idx not in Candidate_Out_Edge_Cluster_List:
				Candidate_Out_Edge_Cluster_List.append(src_taxa_clust_idx)
		if (res == -1):
			Cluster_Info_Dict[dest_taxa_clust_idx]._AddPossibleR1(src_taxa_clust_idx)
			if dest_taxa_clust_idx not in Candidate_Out_Edge_Cluster_List:
				Candidate_Out_Edge_Cluster_List.append(dest_taxa_clust_idx)

	"""
	case E - src cluster size = 1 and dest cluster size > 1
	check R1 relation between this cluster pair
	"""
	if (len(src_cluster_taxa_list) == 1) and (len(dest_cluster_taxa_list) > 1):
		res = CheckCandidateR1R2Reln(src_cluster_taxa_list, dest_cluster_taxa_list)
		if (res == 1):
			Cluster_Info_Dict[src_taxa_clust_idx]._AddPossibleR1(dest_taxa_clust_idx)
			if src_taxa_clust_idx not in Candidate_Out_Edge_Cluster_List:
				Candidate_Out_Edge_Cluster_List.append(src_taxa_clust_idx)


	#"""
	#case D - we check whether R4 relation between the cluster pair can be established
	#provided it does not provide any conflict
	#"""
	#if (CheckTransitiveConflict(src_taxa_label, dest_taxa_label, Reachability_Graph_Mat, RELATION_R4, Output_Text_File) == 0):
		#return 1, RELATION_R4

	return 0, RELATION_R4

#-------------------------------------------------------
#""" 
#this function processes individual queues (score lists) to construct the final supertree 
#"""
#def Proc_Queue(Reachability_Graph_Mat, Output_Text_File):
	#Inp_Queue = Cost_List_Taxa_Pair_Multi_Reln
	#while (0 < len(Inp_Queue)):
		#""" 
		#extract the 1st element of "Inp_Queue" 
		#since it is sorted to have max support score value at the beginning 
		#"""
		#outlist = Heap_Extract_Max(Inp_Queue)
		
		#src_taxa_label = outlist[0]
		#dest_taxa_label = outlist[1]
		#reln_type = outlist[2]
		#conn_score = TaxaPair_Reln_Dict[(src_taxa_label, dest_taxa_label)]._GetEdgeCost_ConnReln(reln_type)
		
		#if (DEBUG_LEVEL >= 2):
			#fp = open(Output_Text_File, 'a')
			#fp.write('\n ===>> CONFLICTING QUEUE -- ')      
			#fp.write(' current extracted max element: ' + str(src_taxa_label) + ' and ' + str(dest_taxa_label) + \
								#' reln type: ' + str(reln_type) + ' conn score: ' + str(conn_score))
			#fp.close()

		##-------------------------------------------------
		### comment - sourya
		
		##""" 
		##if the current support score based relation does not induce a conflict to the existing configuration of the final supertree
		##then include the current relation (and resolve corresponding couplet) in it 
		##"""
		##conflict_detection = Possible_Conflict_Curr_Reln(src_taxa_label, dest_taxa_label, \
			##Reachability_Graph_Mat, reln_type, Output_Text_File, conn_score)
		
		##if (conflict_detection == 0):
			##""" 
			##current element does not create a cycle / conflict
			##not that it is already present in the supertree 
			##valid connection is found - append this connection to the final formed tree 
			##"""
			##if (DEBUG_LEVEL >= 2):
				##fp = open(Output_Text_File, 'a')    
				##queue_str = 'CONFLICTING QUEUE (lower priority)'
				##fp.write('\n ==>>>>>>>>> NEW CONN --- ' + str(queue_str) + 'nodes to be connected: ' + str(src_taxa_label) + ' and ' + str(dest_taxa_label) + \
							##' nodes indices: ' + str(COMPLETE_INPUT_TAXA_LIST.index(src_taxa_label)) + ' and ' + str(COMPLETE_INPUT_TAXA_LIST.index(dest_taxa_label)) + \
							##' relation type: ' + str(reln_type) + ' conn score: ' + str(conn_score))
				##fp.close()
			
			### also update the reachability graph information
			##Reachability_Graph_Mat = AdjustReachGraph(Reachability_Graph_Mat, src_taxa_label, dest_taxa_label, reln_type)
			
		### end comment - sourya
		##-------------------------------------------------
		### add - sourya
		
		##""" 
		##if the current support score based relation does not induce a conflict to 
		##the existing configuration of the final supertree
		##then include the current relation (and resolve corresponding couplet) in it 
		##"""
		##conflict_detection = Possible_Conflict_Curr_Reln(src_taxa_label, dest_taxa_label, \
			##Reachability_Graph_Mat, reln_type, Output_Text_File, conn_score)
		
		##if (conflict_detection == 0):
			##"""
			##Case 1 --- we check whether the support score is positive
			##in such a case, we immediately process the relation
			##"""
			##if (conn_score > 0):
				### process the relation and add it in the list of supported relations
				##if (DEBUG_LEVEL >= 2):
					##fp = open(Output_Text_File, 'a')    
					##queue_str = 'CONFLICTING QUEUE (lower priority)'
					##fp.write('\n ==>>>>>>>>> NEW CONN --- ' + str(queue_str) + 'nodes to be connected: ' \
						##+ str(src_taxa_label) + ' and ' + str(dest_taxa_label) + \
								##' nodes indices: ' + str(COMPLETE_INPUT_TAXA_LIST.index(src_taxa_label)) \
									##+ ' and ' + str(COMPLETE_INPUT_TAXA_LIST.index(dest_taxa_label)) + \
								##' relation type: ' + str(reln_type) + ' conn score: ' + str(conn_score))
					##fp.close()
				
				### also update the reachability graph information
				##Reachability_Graph_Mat = AdjustReachGraph(Reachability_Graph_Mat, src_taxa_label, dest_taxa_label, reln_type)
			##else:
				##"""
				##case 2 - here  the support score is negative
				##so we have to check whether the relation can be really included
				##"""
				##if (CheckConnectionPossible(src_taxa_label, dest_taxa_label, reln_type, Output_Text_File) == True):
					### process the relation and add it in the list of supported relations
					##if (DEBUG_LEVEL >= 2):
						##fp = open(Output_Text_File, 'a')    
						##queue_str = 'CONFLICTING QUEUE (lower priority)'
						##fp.write('\n ==>>>>>>>>> NEW CONN --- ' + str(queue_str) + 'nodes to be connected: ' \
							##+ str(src_taxa_label) + ' and ' + str(dest_taxa_label) + \
									##' nodes indices: ' + str(COMPLETE_INPUT_TAXA_LIST.index(src_taxa_label)) \
										##+ ' and ' + str(COMPLETE_INPUT_TAXA_LIST.index(dest_taxa_label)) + \
									##' relation type: ' + str(reln_type) + ' conn score: ' + str(conn_score))
						##fp.close()

					### also update the reachability graph information
					##Reachability_Graph_Mat = AdjustReachGraph(Reachability_Graph_Mat, src_taxa_label, dest_taxa_label, reln_type)
				
				##else:
					##""" 
					##here the relation can not be added
					##so remove this relation from the set of allowed relation for this couplet
					##"""
					##l = (src_taxa_label, dest_taxa_label)
					##TaxaPair_Reln_Dict[l]._RemoveAllowedReln(reln_type)
				
		##else:
			##"""
			##current relation produces a conflict
			##so remove this relation from the set of allowed relation for this couplet
			##"""
			##l = (src_taxa_label, dest_taxa_label)
			##TaxaPair_Reln_Dict[l]._RemoveAllowedReln(reln_type)

		### end add - sourya
		##-------------------------------------------------
		### add - sourya
		
		##""" 
		##if the couplet is already connected, as shown by the entries in Reachability_Graph_Mat
		##then there is no need for any connection
		##"""
		##existing_conn = CheckExistingConn(src_taxa_label, dest_taxa_label, Reachability_Graph_Mat, reln_type, Output_Text_File)
		##if (existing_conn == 0):
			##"""
			##there is no apparent existing relationship between the couplet
			##"""
			##"""
			##Case 1 --- we check whether the support score is positive
			##in such a case, we immediately process the relation
			##provided the relation is not conflicting
			##"""
			##if (conn_score > 0):
				##if (CheckTransitiveConflict(src_taxa_label, dest_taxa_label, Reachability_Graph_Mat, reln_type, Output_Text_File) == 0):
					### process the relation and add it in the list of supported relations
					##if (DEBUG_LEVEL >= 2):
						##fp = open(Output_Text_File, 'a')    
						##queue_str = 'CONFLICTING QUEUE (lower priority)'
						##fp.write('\n ==>>>>>>>>> NEW CONN --- ' + str(queue_str) + 'nodes to be connected: ' \
							##+ str(src_taxa_label) + ' and ' + str(dest_taxa_label) + \
									##' nodes indices: ' + str(COMPLETE_INPUT_TAXA_LIST.index(src_taxa_label)) \
										##+ ' and ' + str(COMPLETE_INPUT_TAXA_LIST.index(dest_taxa_label)) + \
									##' relation type: ' + str(reln_type) + ' conn score: ' + str(conn_score))
						##fp.close()
					
					### also update the reachability graph information
					##Reachability_Graph_Mat = AdjustReachGraph(Reachability_Graph_Mat, src_taxa_label, dest_taxa_label, reln_type)

			##else:
				##"""
				##Case 2 --- here the support score is negative
				##in such a case, we check the relation which best suits the couplet
				##provided the relation is not conflicting
				##"""
				##conn_done, connecting_reln_type = ReturnConnectingReln(src_taxa_label, dest_taxa_label, Reachability_Graph_Mat, Output_Text_File)
				##if (conn_done == 1):
					##"""
					##the couplet is now related with the 'connecting_reln_type'
					##update the Reachability_Graph_Mat also
					##"""
					##if (DEBUG_LEVEL >= 2):
						##fp = open(Output_Text_File, 'a')    
						##queue_str = 'CONFLICTING QUEUE (lower priority)'
						##fp.write('\n ==>>>>>>>>> NEW CONN --- ' + str(queue_str) + 'nodes to be connected: ' \
							##+ str(src_taxa_label) + ' and ' + str(dest_taxa_label) + \
									##' nodes indices: ' + str(COMPLETE_INPUT_TAXA_LIST.index(src_taxa_label)) \
										##+ ' and ' + str(COMPLETE_INPUT_TAXA_LIST.index(dest_taxa_label)) + \
									##' relation type: ' + str(connecting_reln_type))
						##fp.close()
					
					### also update the reachability graph information
					##Reachability_Graph_Mat = AdjustReachGraph(Reachability_Graph_Mat, src_taxa_label, dest_taxa_label, connecting_reln_type)

				##else:
					##"""
					##no relation can be established between the couplet
					##"""
					##if (DEBUG_LEVEL >= 2):
						##fp = open(Output_Text_File, 'a')    
						##fp.write('\n Although the couplet is not connected, no relation can be established between them')
						##fp.close()
		
		##else:
			##"""
			##already the couplet is connected - so print that information
			##"""
			##if (DEBUG_LEVEL >= 2):
				##fp = open(Output_Text_File, 'a')    
				##fp.write('\n The couplet (' + str(src_taxa_label) + ', ' + str(dest_taxa_label) + ') is already connected') 
				##fp.close()
			
		### end add - sourya
		##-------------------------------------------------
		
		## this was used at the last code - sourya
		
		### add - sourya
		##""" 
		##if the couplet is already connected, as shown by the entries in Reachability_Graph_Mat
		##then there is no need for any connection
		##"""
		##existing_conn = CheckExistingConn(src_taxa_label, dest_taxa_label, Reachability_Graph_Mat, reln_type, Output_Text_File)
		##if (existing_conn == 0):
			##"""
			##there is no apparent existing relationship between the couplet
			##"""
			##"""
			##Case 1 --- we check whether the support score is positive
			##in such a case, we immediately process the relation
			##provided the relation is not conflicting
			##"""
			##if (conn_score > 0):
				##if (CheckTransitiveConflict(src_taxa_label, dest_taxa_label, Reachability_Graph_Mat, reln_type, Output_Text_File) == 0):
					### process the relation and add it in the list of supported relations
					##if (DEBUG_LEVEL >= 2):
						##fp = open(Output_Text_File, 'a')    
						##queue_str = 'CONFLICTING QUEUE (lower priority)'
						##fp.write('\n ==>>>>>>>>> NEW CONN --- ' + str(queue_str) + 'nodes to be connected: ' \
							##+ str(src_taxa_label) + ' and ' + str(dest_taxa_label) + \
									##' nodes indices: ' + str(COMPLETE_INPUT_TAXA_LIST.index(src_taxa_label)) \
										##+ ' and ' + str(COMPLETE_INPUT_TAXA_LIST.index(dest_taxa_label)) + \
									##' relation type: ' + str(reln_type) + ' conn score: ' + str(conn_score))
						##fp.close()
					
					### also update the reachability graph information
					##Reachability_Graph_Mat = AdjustReachGraph(Reachability_Graph_Mat, src_taxa_label, dest_taxa_label, reln_type)
					
			##else:
				##"""
				##Case 2 --- here the support score is negative
				##in such a case, we first check if R3 can be considered as a relation between them
				##otherwise, proceed with the current relation
				##"""
				##l = (src_taxa_label, dest_taxa_label)
				##reln_list = TaxaPair_Reln_Dict[l]._GetAllowedRelnList()
				##conn_done = 0	# serves as the flag variable
				
				##if RELATION_R3 in reln_list:
					##if (TaxaPair_Reln_Dict[l]._Check_Reln_R3_Majority(Output_Text_File) == True):
						##if (CheckTransitiveConflict(src_taxa_label, dest_taxa_label, Reachability_Graph_Mat, RELATION_R3, Output_Text_File) == 0):
							##"""
							##the couplet is now related with the 'RELATION_R3'
							##update the Reachability_Graph_Mat also
							##"""
							##if (DEBUG_LEVEL >= 2):
								##fp = open(Output_Text_File, 'a')    
								##queue_str = 'CONFLICTING QUEUE (lower priority)'
								##fp.write('\n ==>>>>>>>>> NEW CONN --- ' + str(queue_str) + 'nodes to be connected: ' \
									##+ str(src_taxa_label) + ' and ' + str(dest_taxa_label) + \
											##' nodes indices: ' + str(COMPLETE_INPUT_TAXA_LIST.index(src_taxa_label)) \
												##+ ' and ' + str(COMPLETE_INPUT_TAXA_LIST.index(dest_taxa_label)) + \
											##' relation type: ' + str(RELATION_R3))
								##fp.close()
							
							### also update the reachability graph information
							##Reachability_Graph_Mat = AdjustReachGraph(Reachability_Graph_Mat, src_taxa_label, dest_taxa_label, RELATION_R3)
							
							##conn_done = 1	# update the flag variable
				
				##"""
				##otherwise, we check whether the current relation can be included 
				##"""
				##if (conn_done == 0):
					##if (CheckTransitiveConflict(src_taxa_label, dest_taxa_label, Reachability_Graph_Mat, reln_type, Output_Text_File) == 0):
						### process the relation and add it in the list of supported relations
						##if (DEBUG_LEVEL >= 2):
							##fp = open(Output_Text_File, 'a')    
							##queue_str = 'CONFLICTING QUEUE (lower priority)'
							##fp.write('\n ==>>>>>>>>> NEW CONN --- ' + str(queue_str) + 'nodes to be connected: ' \
								##+ str(src_taxa_label) + ' and ' + str(dest_taxa_label) + \
										##' nodes indices: ' + str(COMPLETE_INPUT_TAXA_LIST.index(src_taxa_label)) \
											##+ ' and ' + str(COMPLETE_INPUT_TAXA_LIST.index(dest_taxa_label)) + \
										##' relation type: ' + str(reln_type) + ' conn score: ' + str(conn_score))
							##fp.close()
						
						### also update the reachability graph information
						##Reachability_Graph_Mat = AdjustReachGraph(Reachability_Graph_Mat, src_taxa_label, dest_taxa_label, reln_type)
		##else:
			##"""
			##already the couplet is connected - so print that information
			##"""
			##if (DEBUG_LEVEL >= 2):
				##fp = open(Output_Text_File, 'a')    
				##fp.write('\n The couplet (' + str(src_taxa_label) + ', ' + str(dest_taxa_label) + ') is already connected') 
				##fp.close()
			
		### end add - sourya
		
		##-------------------------------------------------
		## new code - sourya - 16.02.2016

		#"""
		#sourya - we modify this code
		#we do not insert the R3 relation in the queue
		#"""
		#if (reln_type == RELATION_R3):
			#continue

		#""" 
		#if the couplet is already connected, as shown by the entries in Reachability_Graph_Mat
		#then there is no need for any connection
		#"""
		
		#existing_conn = CheckExistingConn(src_taxa_label, dest_taxa_label, Reachability_Graph_Mat, reln_type, Output_Text_File)
		#if (existing_conn == 0):
			#"""
			#there is no apparent existing relationship between the couplet
			#"""
			#"""
			#Case 1A - relation is either R1 or R2
			#in such a case, immediately apply the relation
			#provided the relation is not conflicting
			#"""
			#if (reln_type == RELATION_R1) or (reln_type == RELATION_R2):
				#"""
				#first, check whether the relation induces conflict
				#otherwise, we 
				#"""
				#if (CheckTransitiveConflict(src_taxa_label, dest_taxa_label, Reachability_Graph_Mat, reln_type, Output_Text_File) == 0):
					#if (DEBUG_LEVEL >= 2):
						#fp = open(Output_Text_File, 'a')    
						#fp.write('\n ==>>>>>>>>> NEW CONN --- nodes to be connected: ' \
							#+ str(src_taxa_label) + ' and ' + str(dest_taxa_label) + \
									#' nodes indices: ' + str(COMPLETE_INPUT_TAXA_LIST.index(src_taxa_label)) \
										#+ ' and ' + str(COMPLETE_INPUT_TAXA_LIST.index(dest_taxa_label)) + \
									#' relation type: ' + str(reln_type) + ' conn score: ' + str(conn_score))
						#fp.close()
					
					## also update the reachability graph information
					#Reachability_Graph_Mat = AdjustReachGraph(Reachability_Graph_Mat, src_taxa_label, dest_taxa_label, reln_type)
					
			#"""
			#case 1B - relation is R4
			#"""
			#if (reln_type == RELATION_R4):
				#"""
				#here we check about possible R1 or R2 relation between the cluster pairs
				#otherwise, we use the relation R4
				#"""
				#conn_done, target_reln = CheckHiddenR1R2Reln(Reachability_Graph_Mat, src_taxa_label, dest_taxa_label, Output_Text_File)
				#if (conn_done == 1):
					#if (DEBUG_LEVEL >= 2):
						#fp = open(Output_Text_File, 'a')    
						#fp.write('\n ==>>>>>>>>> NEW CONN --- nodes to be connected: ' \
							#+ str(src_taxa_label) + ' and ' + str(dest_taxa_label) + \
							#' nodes indices: ' + str(COMPLETE_INPUT_TAXA_LIST.index(src_taxa_label)) \
							#+ ' and ' + str(COMPLETE_INPUT_TAXA_LIST.index(dest_taxa_label)) + \
							#' relation type: ' + str(target_reln) + ' conn score: ' \
							#+ str(TaxaPair_Reln_Dict[(src_taxa_label, dest_taxa_label)]._GetEdgeCost_ConnReln(target_reln)))
						#fp.close()

				## also update the reachability graph information
				#Reachability_Graph_Mat = AdjustReachGraph(Reachability_Graph_Mat, src_taxa_label, dest_taxa_label, target_reln)
							
		#else:
			#"""
			#already the couplet is connected - so print that information
			#"""
			#if (DEBUG_LEVEL >= 2):
				#fp = open(Output_Text_File, 'a')    
				#fp.write('\n The couplet (' + str(src_taxa_label) + ', ' + str(dest_taxa_label) + ') is already connected') 
				#fp.close()
					
		## end new code - sourya - 16.02.2016 
		##-------------------------------------------------
		
	#return Reachability_Graph_Mat
  