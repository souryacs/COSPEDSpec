#!/usr/bin/env python

import Header
from Header import *
import Cost_Update
from Cost_Update import *
import ReachGraph_Update
from ReachGraph_Update import *
import Conflict_Detect
from Conflict_Detect import *

#---------------------------------------------
"""
this function establishes the R1 relation from the src cluster to the dest cluster
"""
def EstablishR1Reln(src_clust, dest_clust, Reachability_Graph_Mat, Output_Text_File):
	if (DEBUG_LEVEL >= 2):
		fp = open(Output_Text_File, 'a')
		fp.write('\n New directed R1 edge from the cluster ' + str(src_clust) + '  to the cluster : ' + str(dest_clust))
		fp.close()
	# update the reachability graph
	Reachability_Graph_Mat = AdjustReachGraph(Reachability_Graph_Mat, src_clust, dest_clust, RELATION_R1)
	return Reachability_Graph_Mat

#---------------------------------------------
"""
this function processes all clusters having candidate out edge information
"""
def Process_Candidate_Out_Edge_Cluster_List(Reachability_Graph_Mat, DIST_MAT_TYPE, Output_Text_File):
	for cl in Candidate_Out_Edge_Cluster_List:
		if (DEBUG_LEVEL >= 2):
			fp = open(Output_Text_File, 'a')
			fp.write('\n ********** current cl no: ' + str(cl)) 
			fp.close()

		"""
		check if the cluster cl has any directed out edge already
		otherwise, place all its candidate R1 clusters as a child to its parent node
		"""
		if (Cluster_Info_Dict[cl]._Get_Outdegree() == 0):
			"""
			here, assign out edge from all the parent clusters of cl to this cluster Y
			"""
			for parent_cl in Cluster_Info_Dict[cl]._GetInEdgeList():
				for Y in Cluster_Info_Dict[cl]._GetPossibleR1List():
					if (CheckExistingConn(parent_cl, Y, Reachability_Graph_Mat, RELATION_R1, Output_Text_File) == 0):
						if (CheckTransitiveConflict(parent_cl, Y, Reachability_Graph_Mat, RELATION_R1, Output_Text_File) == 0):
							Reachability_Graph_Mat = EstablishR1Reln(parent_cl, Y, Reachability_Graph_Mat, Output_Text_File)
		else:
		
			"""
			create a list with individual elements having two fields
			one is the cluster index
			another is the XL value of the cluster with cl
			"""
			temp_subl = []
			for Y in Cluster_Info_Dict[cl]._GetPossibleR1List():
				"""
				excess gene count between the cluster Y and the cluster cl
				"""
				xl_cl_Y = FindAvgXL(Cluster_Info_Dict[cl]._GetSpeciesList(),\
					Cluster_Info_Dict[Y]._GetSpeciesList(), DIST_MAT_TYPE, 1)
				
				elem = [Y, xl_cl_Y]
				temp_subl.append(elem)
				
			"""
			sort the temp_subl with the second field (XL) as key element
			in ascending order
			"""
			temp_subl.sort(key=lambda x: x[1])
		
			"""
			for each Y belonging to the temp_subl (already sorted)
			1) compute XL(cl, Y)
			explore individual children "child_cl" of the cluster cl
			compute XL(Y, child_cl)
			and also compute XL(cl, child_cl)
			"""
			for Y_idx in range(len(temp_subl)):
				Y = temp_subl[Y_idx][0]
				xl_cl_Y = temp_subl[Y_idx][1]
				
				"""
				analyze Y only if there is no already existing connection between cl and Y
				"""
				if (CheckExistingConn(cl, Y, Reachability_Graph_Mat, RELATION_R1, Output_Text_File) == 0):
					if (CheckTransitiveConflict(cl, Y, Reachability_Graph_Mat, RELATION_R1, Output_Text_File) == 0):
						if (DEBUG_LEVEL >= 2):
							fp = open(Output_Text_File, 'a')
							fp.write('\n Processing possible R1 cluster : ' + str(Y)) 
							fp.close()
						
						if (DEBUG_LEVEL >= 2):
							fp = open(Output_Text_File, 'a')
							fp.write('\n FINAL Average excess gene count between ' + str(cl) \
								+ '  and the cluster ' + str(Y) + ' is: ' + str(xl_cl_Y)) 
							fp.close()
						
						"""
						this is the average of excess gene count between Y and every child of cl
						"""
						xl_Y_childcl = 0
						"""
						this is the average of excess gene count between cl and every child of cl
						"""
						xl_cl_childcl = 0
						
						for child_cl in Cluster_Info_Dict[cl]._GetOutEdgeList():
							curr_xl_Y_childcl = FindAvgXL(Cluster_Info_Dict[child_cl]._GetSpeciesList(), \
								Cluster_Info_Dict[Y]._GetSpeciesList(), DIST_MAT_TYPE, 1)
							xl_Y_childcl = xl_Y_childcl + curr_xl_Y_childcl
							
							curr_xl_cl_childcl = FindAvgXL(Cluster_Info_Dict[child_cl]._GetSpeciesList(), \
								Cluster_Info_Dict[cl]._GetSpeciesList(), DIST_MAT_TYPE, 1)
							xl_cl_childcl = xl_cl_childcl + curr_xl_cl_childcl
							
							if (DEBUG_LEVEL >= 2):
								fp = open(Output_Text_File, 'a')
								fp.write('\n excess gene count between (child) ' + str(child_cl) \
									+ '  and the cluster ' + str(Y) + ' is: ' + str(curr_xl_Y_childcl)) 
								fp.write('\n excess gene count between (child) ' + str(child_cl) \
									+ '  and the cluster ' + str(cl) + ' is: ' + str(curr_xl_cl_childcl)) 
								fp.close()
						
						"""
						we average the XL measures with respect to all the clusters
						"""
						xl_Y_childcl = (xl_Y_childcl * 1.0) / len(Cluster_Info_Dict[cl]._GetOutEdgeList())
						xl_cl_childcl = (xl_cl_childcl * 1.0) / len(Cluster_Info_Dict[cl]._GetOutEdgeList())

						if (DEBUG_LEVEL >= 2):
							fp = open(Output_Text_File, 'a')
							fp.write('\n FINAL Avg excess gene count between (child) clusters and the cluster ' \
								+ str(Y) + ' is: ' + str(xl_Y_childcl)) 
							fp.write('\n FINAL Avg excess gene count between clusters and the cluster ' \
								+ str(cl) + ' is: ' + str(xl_cl_childcl)) 
							fp.close()

						"""
						the condition for topology (A,(B,C)) is that XL(B,C) should be less than both XL(A,B) and XL(A,C)
						"""
						if (xl_Y_childcl < xl_cl_Y) and (xl_Y_childcl < xl_cl_childcl):
							"""
							Y can be placed as the child of cl
							"""
							Reachability_Graph_Mat = EstablishR1Reln(cl, Y, Reachability_Graph_Mat, Output_Text_File)
						else:
							"""
							here, assign out edge from all the parent clusters of cl to this cluster Y
							"""
							for parent_cl in Cluster_Info_Dict[cl]._GetInEdgeList():
								if (CheckExistingConn(parent_cl, Y, Reachability_Graph_Mat, RELATION_R1, Output_Text_File) == 0):
									if (CheckTransitiveConflict(parent_cl, Y, Reachability_Graph_Mat, RELATION_R1, Output_Text_File) == 0):
										Reachability_Graph_Mat = EstablishR1Reln(parent_cl, Y, Reachability_Graph_Mat, Output_Text_File)
					
	return Reachability_Graph_Mat


##---------------------------------------------
#"""
#this function processes all clusters having candidate out edge information
#"""
#def Process_Candidate_Out_Edge_Cluster_List(Reachability_Graph_Mat, DIST_MAT_TYPE, Output_Text_File):
	#for cl in Candidate_Out_Edge_Cluster_List:
		#if (DEBUG_LEVEL >= 2):
			#fp = open(Output_Text_File, 'a')
			#fp.write('\n ********** current cl no: ' + str(cl)) 
			#fp.close()
		
		#"""
		#check if the cluster cl has any directed out edge already
		#otherwise, place all its candidate R1 clusters as a child to its parent node
		#"""
		#if (Cluster_Info_Dict[cl]._Get_Outdegree() == 0):
			
			### comment - sourya
			##"""
			##for each cluster Y in the possible candidate R1 cluster list of the cluster cl
			##check the in edge list 'In' of Y
			##construct a list "Y_In" which contains all the clusters p in 'In'
			##such that p and cl are not related by any means
			##if the size of list "Y_In" is > 0, compute minXL = XL(p, Y) for any p in "Y_In"
			##if minXL > XL(cl, Y) then add Y as a descendant of cl
			##else (or even if len("Y_In") == 0), add Y as a child to individual parent clusters of cl
			##"""
			##for Y in Cluster_Info_Dict[cl]._GetPossibleR1List():
				##"""
				##analyze Y only if there is no already existing connection between cl and Y
				##"""
				##if (CheckExistingConn(cl, Y, Reachability_Graph_Mat, RELATION_R1, Output_Text_File) == 0):
					##if (CheckTransitiveConflict(cl, Y, Reachability_Graph_Mat, RELATION_R1, Output_Text_File) == 0):
						##if (DEBUG_LEVEL >= 2):
							##fp = open(Output_Text_File, 'a')
							##fp.write('\n Processing possible R1 cluster : ' + str(Y)) 
							##fp.close()
						##"""
						##construct "Y_In", as mentioned above
						##also construct a list to contain XL(p, Y) values
						##"""
						##Y_In_Clust_List = []
						##XL_Y_In_Clust_List = []
						##for p in Cluster_Info_Dict[Y]._GetInEdgeList():
							##if (CheckExistingConn(cl, p, Reachability_Graph_Mat, RELATION_R1, Output_Text_File) == 0):
								##"""
								##update "Y_In" and the XL list
								##"""
								##Y_In_Clust_List.append(p)
								##XL_Y_In_Clust_List.append(FindAvgXL(Cluster_Info_Dict[p]._GetSpeciesList(), \
									##Cluster_Info_Dict[Y]._GetSpeciesList(), DIST_MAT_TYPE, 1))
							
								##if (len(Y_In_Clust_List) > 0):
									##"""
									##minimum of XL values within min_XL_Y_In_Clust_List
									##"""
									##min_XL_Y_In_Clust_List = min(XL_Y_In_Clust_List)
									##"""
									##excess gene count between the cluster Y and the cluster cl
									##"""
									##xl_cl_Y = FindAvgXL(Cluster_Info_Dict[cl]._GetSpeciesList(), \
										##Cluster_Info_Dict[Y]._GetSpeciesList(), DIST_MAT_TYPE, 1)
									
									##if (xl_cl_Y < min_XL_Y_In_Clust_List):
										### add Y as a descendant of cl
										##Reachability_Graph_Mat = EstablishR1Reln(cl, Y, Reachability_Graph_Mat, Output_Text_File)
									##else:
										##"""
										##here, assign out edge from all the parent clusters of cl to this cluster Y
										##"""
										##for parent_cl in Cluster_Info_Dict[cl]._GetInEdgeList():
											##if (CheckExistingConn(parent_cl, Y, Reachability_Graph_Mat, RELATION_R1, Output_Text_File) == 0):
												##if (CheckTransitiveConflict(parent_cl, Y, Reachability_Graph_Mat, RELATION_R1, Output_Text_File) == 0):
													##Reachability_Graph_Mat = EstablishR1Reln(parent_cl, Y, Reachability_Graph_Mat, Output_Text_File)
								##else:
									##"""
									##here, assign out edge from all the parent clusters of cl to this cluster Y
									##"""
									##for parent_cl in Cluster_Info_Dict[cl]._GetInEdgeList():
										##if (CheckExistingConn(parent_cl, Y, Reachability_Graph_Mat, RELATION_R1, Output_Text_File) == 0):
											##if (CheckTransitiveConflict(parent_cl, Y, Reachability_Graph_Mat, RELATION_R1, Output_Text_File) == 0):
												##Reachability_Graph_Mat = EstablishR1Reln(parent_cl, Y, Reachability_Graph_Mat, Output_Text_File)
			
			### end comment - sourya

			#"""
			#here, assign out edge from all the parent clusters of cl to this cluster Y
			#"""
			#for parent_cl in Cluster_Info_Dict[cl]._GetInEdgeList():
				#for Y in Cluster_Info_Dict[cl]._GetPossibleR1List():
					#if (CheckExistingConn(parent_cl, Y, Reachability_Graph_Mat, RELATION_R1, Output_Text_File) == 0):
						#if (CheckTransitiveConflict(parent_cl, Y, Reachability_Graph_Mat, RELATION_R1, Output_Text_File) == 0):
							#Reachability_Graph_Mat = EstablishR1Reln(parent_cl, Y, Reachability_Graph_Mat, Output_Text_File)
			
		#else:
			#"""
			#for each Y belonging to the candidate R1 list of the cluster cl
			#1) compute XL(cl, Y)
			#explore individual children "child_cl" of the cluster cl
			#compute XL(Y, child_cl)
			#and also compute XL(cl, child_cl)
			#"""
			#for Y in Cluster_Info_Dict[cl]._GetPossibleR1List():
				#"""
				#analyze Y only if there is no already existing connection between cl and Y
				#"""
				#if (CheckExistingConn(cl, Y, Reachability_Graph_Mat, RELATION_R1, Output_Text_File) == 0):
					#if (CheckTransitiveConflict(cl, Y, Reachability_Graph_Mat, RELATION_R1, Output_Text_File) == 0):

						#if (DEBUG_LEVEL >= 2):
							#fp = open(Output_Text_File, 'a')
							#fp.write('\n Processing possible R1 cluster : ' + str(Y)) 
							#fp.close()
						
						#"""
						#excess gene count between the cluster Y and the cluster cl
						#"""
						#xl_cl_Y = FindAvgXL(Cluster_Info_Dict[cl]._GetSpeciesList(),\
							#Cluster_Info_Dict[Y]._GetSpeciesList(), DIST_MAT_TYPE, 1)
						
						#if (DEBUG_LEVEL >= 2):
							#fp = open(Output_Text_File, 'a')
							#fp.write('\n FINAL Average excess gene count between ' + str(cl) \
								#+ '  and the cluster ' + str(Y) + ' is: ' + str(xl_cl_Y)) 
							#fp.close()
						
						#"""
						#this is the average of excess gene count between Y and every child of cl
						#"""
						#xl_Y_childcl = 0
						#"""
						#this is the average of excess gene count between cl and every child of cl
						#"""
						#xl_cl_childcl = 0
						
						#for child_cl in Cluster_Info_Dict[cl]._GetOutEdgeList():
							#curr_xl_Y_childcl = FindAvgXL(Cluster_Info_Dict[child_cl]._GetSpeciesList(), \
								#Cluster_Info_Dict[Y]._GetSpeciesList(), DIST_MAT_TYPE, 1)
							#xl_Y_childcl = xl_Y_childcl + curr_xl_Y_childcl
							
							#curr_xl_cl_childcl = FindAvgXL(Cluster_Info_Dict[child_cl]._GetSpeciesList(), \
								#Cluster_Info_Dict[cl]._GetSpeciesList(), DIST_MAT_TYPE, 1)
							#xl_cl_childcl = xl_cl_childcl + curr_xl_cl_childcl
							
							#if (DEBUG_LEVEL >= 2):
								#fp = open(Output_Text_File, 'a')
								#fp.write('\n excess gene count between (child) ' + str(child_cl) \
									#+ '  and the cluster ' + str(Y) + ' is: ' + str(curr_xl_Y_childcl)) 
								#fp.write('\n excess gene count between (child) ' + str(child_cl) \
									#+ '  and the cluster ' + str(cl) + ' is: ' + str(curr_xl_cl_childcl)) 
								#fp.close()
						
						#"""
						#we average the XL measures with respect to all the clusters
						#"""
						#xl_Y_childcl = (xl_Y_childcl * 1.0) / len(Cluster_Info_Dict[cl]._GetOutEdgeList())
						#xl_cl_childcl = (xl_cl_childcl * 1.0) / len(Cluster_Info_Dict[cl]._GetOutEdgeList())

						#if (DEBUG_LEVEL >= 2):
							#fp = open(Output_Text_File, 'a')
							#fp.write('\n FINAL Avg excess gene count between (child) clusters and the cluster ' \
								#+ str(Y) + ' is: ' + str(xl_Y_childcl)) 
							#fp.write('\n FINAL Avg excess gene count between clusters and the cluster ' \
								#+ str(cl) + ' is: ' + str(xl_cl_childcl)) 
							#fp.close()

						#"""
						#the condition for topology (A,(B,C)) is that XL(B,C) should be less than both XL(A,B) and XL(A,C)
						#"""
						#if (xl_Y_childcl < xl_cl_Y) and (xl_Y_childcl < xl_cl_childcl):
							#"""
							#Y can be placed as the child of cl
							#"""
							#Reachability_Graph_Mat = EstablishR1Reln(cl, Y, Reachability_Graph_Mat, Output_Text_File)
						#else:
							#"""
							#here, assign out edge from all the parent clusters of cl to this cluster Y
							#"""
							#for parent_cl in Cluster_Info_Dict[cl]._GetInEdgeList():
								#if (CheckExistingConn(parent_cl, Y, Reachability_Graph_Mat, RELATION_R1, Output_Text_File) == 0):
									#if (CheckTransitiveConflict(parent_cl, Y, Reachability_Graph_Mat, RELATION_R1, Output_Text_File) == 0):
										#Reachability_Graph_Mat = EstablishR1Reln(parent_cl, Y, Reachability_Graph_Mat, Output_Text_File)
					
	#return Reachability_Graph_Mat

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
		src_taxa_clust_idx = Taxa_Info_Dict[src_taxa_label]._Get_Taxa_Part_Clust_Idx()
		dest_taxa_clust_idx = Taxa_Info_Dict[dest_taxa_label]._Get_Taxa_Part_Clust_Idx()
		
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
				if (CheckExistingConn(src_taxa_clust_idx, dest_taxa_clust_idx, \
							Reachability_Graph_Mat, reln_type, Output_Text_File) == 0):
					"""
					there is no apparent existing relationship between the couplet
					"""
					"""
					if the R1 / R2 relation is strictly consensus then we apply the relation directly
					otherwise, we insert the relation in the candidate set of relations
					"""
					if (conn_score <= 0):
						conn_done, target_reln = CheckConsensusR1R2NegScoreCase(reln_type, Reachability_Graph_Mat, \
							src_taxa_clust_idx, dest_taxa_clust_idx, Output_Text_File)
						#conn_done, target_reln = CheckHiddenR1R2Reln(Reachability_Graph_Mat, src_taxa_clust_idx, dest_taxa_clust_idx, Output_Text_File)
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

							"""
							# also update the reachability graph information
							Note: This updation only takes place if conn_done = 1
							"""
							Reachability_Graph_Mat = AdjustReachGraph(Reachability_Graph_Mat, src_taxa_clust_idx, dest_taxa_clust_idx, target_reln)
					
					else:
						"""
						first, check whether the relation induces conflict
						otherwise, we apply the relations
						"""
						if (CheckTransitiveConflict(src_taxa_clust_idx, dest_taxa_clust_idx, \
								Reachability_Graph_Mat, reln_type, Output_Text_File) == 0):
							if (DEBUG_LEVEL >= 2):
								fp = open(Output_Text_File, 'a')    
								fp.write('\n ==>>>>>>>>> NEW CONN --- nodes to be connected: ' \
									+ str(src_taxa_label) + ' and ' + str(dest_taxa_label) + \
											' nodes indices: ' + str(COMPLETE_INPUT_TAXA_LIST.index(src_taxa_label)) \
												+ ' and ' + str(COMPLETE_INPUT_TAXA_LIST.index(dest_taxa_label)) + \
											' relation type: ' + str(reln_type) + ' conn score: ' + str(conn_score))
								fp.close()
							
							# also update the reachability graph information
							Reachability_Graph_Mat = AdjustReachGraph(Reachability_Graph_Mat, src_taxa_clust_idx, dest_taxa_clust_idx, reln_type)
				
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
				if (CheckExistingConn(src_taxa_clust_idx, dest_taxa_clust_idx, \
							Reachability_Graph_Mat, reln_type, Output_Text_File) == 0):
					"""
					there is no apparent existing relationship between the couplet
					"""
					"""
					here we check about possible R1 or R2 relation between the cluster pairs
					otherwise we do not process the couplet for now
					"""
					conn_done, target_reln = CheckHiddenR1R2Reln(Reachability_Graph_Mat, src_taxa_clust_idx, dest_taxa_clust_idx, Output_Text_File)
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

						"""
						# also update the reachability graph information
						Note: This updation only takes place if conn_done = 1
						"""
						Reachability_Graph_Mat = AdjustReachGraph(Reachability_Graph_Mat, src_taxa_clust_idx, dest_taxa_clust_idx, target_reln)

	return Reachability_Graph_Mat
#-------------------------------------------------------
"""
this function checks for the couplets having R1 / R2 relation as consensus
but with negative support score
"""
def CheckConsensusR1R2NegScoreCase(target_reln_type, Reachability_Graph_Mat, \
	src_taxa_clust_idx, dest_taxa_clust_idx, Output_Text_File):
	"""
	find the taxa list of respective clusters
	"""
	src_cluster_taxa_list = Cluster_Info_Dict[src_taxa_clust_idx]._GetSpeciesList()
	dest_cluster_taxa_list = Cluster_Info_Dict[dest_taxa_clust_idx]._GetSpeciesList()
	if (len(src_cluster_taxa_list) == 1) and (len(dest_cluster_taxa_list) == 1):
		"""
		here the couplet has consensus R1 / R2 relation with respect to the input trees
		"""
		if (target_reln_type == RELATION_R1):
			""" 
			consensus R1 relation from src_taxa_clust_idx to dest_taxa_clust_idx
			"""
			res = Check_Consensus_R1Reln(src_cluster_taxa_list, dest_cluster_taxa_list)
			if (DEBUG_LEVEL >= 2):
				fp = open(Output_Text_File, 'a')
				fp.write('\n ---- src_cluster_taxa_list: ' + str(src_cluster_taxa_list) + ' dest_cluster_taxa_list: ' \
					+ str(dest_cluster_taxa_list) + '  function Check_Consensus_R1Reln -- res: ' + str(res))
				fp.close()
			if (res == 1):
				if (CheckTransitiveConflict(src_taxa_clust_idx, dest_taxa_clust_idx, \
								Reachability_Graph_Mat, RELATION_R1, Output_Text_File) == 0):
					return 1, RELATION_R1

			"""
			otherwise, we insert dest_taxa_clust_idx onto the candidate out edge list of src_taxa_clust_idx
			"""
			Cluster_Info_Dict[src_taxa_clust_idx]._AddPossibleR1(dest_taxa_clust_idx)
			if src_taxa_clust_idx not in Candidate_Out_Edge_Cluster_List:
				Candidate_Out_Edge_Cluster_List.append(src_taxa_clust_idx)
			return 0, RELATION_R4
		
		else:
			""" 
			consensus R2 relation from src_taxa_clust_idx to dest_taxa_clust_idx
			"""
			res = Check_Consensus_R1Reln(dest_cluster_taxa_list, src_cluster_taxa_list)
			if (DEBUG_LEVEL >= 2):
				fp = open(Output_Text_File, 'a')
				fp.write('\n ---- dest_cluster_taxa_list: ' + str(dest_cluster_taxa_list) + ' src_cluster_taxa_list: ' \
					+ str(src_cluster_taxa_list) + '  function Check_Consensus_R1Reln -- res: ' + str(res))
				fp.close()
			if (res == 1):
				if (CheckTransitiveConflict(src_taxa_clust_idx, dest_taxa_clust_idx, \
								Reachability_Graph_Mat, RELATION_R2, Output_Text_File) == 0):
					return 1, RELATION_R2
		
			"""
			otherwise, we insert src_taxa_clust_idx onto the candidate out edge list of dest_taxa_clust_idx
			"""
			Cluster_Info_Dict[dest_taxa_clust_idx]._AddPossibleR1(src_taxa_clust_idx)
			if dest_taxa_clust_idx not in Candidate_Out_Edge_Cluster_List:
				Candidate_Out_Edge_Cluster_List.append(dest_taxa_clust_idx)
			return 0, RELATION_R4

	else:
		"""
		at least one of the clusters have cardinality > 1
		so, we process the default function
		"""
		conn_done, target_reln = CheckHiddenR1R2Reln(Reachability_Graph_Mat, src_taxa_clust_idx, dest_taxa_clust_idx, Output_Text_File)
	
	return conn_done, target_reln

#-------------------------------------------------------
"""
this function checks whether R1 relation from clust1_taxa_list to clust2_taxa_list 
can be established or not, when the consensus relation is R1
This function is called only when both clust1_spec_list and clust2_spec_list 
have cardinality 1, i.e. leaf to leaf R1/R2 relation is considered
"""
def Check_Consensus_R1Reln(clust1_spec_list, clust2_spec_list):

	# explore all taxa pairs of the cluster pair
	for t1 in clust1_spec_list:
		for t2 in clust2_spec_list:
			key1 = (t1, t2)
			key2 = (t2, t1)
			if key1 in TaxaPair_Reln_Dict:
				r1_level_val_ratio = TaxaPair_Reln_Dict[key1]._GetLevelValRatio(0)
				allowed_reln_list = TaxaPair_Reln_Dict[key1]._GetAllowedRelnList()
				r1_freq = TaxaPair_Reln_Dict[key1]._GetEdgeWeight(RELATION_R1)
				r4_freq = TaxaPair_Reln_Dict[key1]._GetEdgeWeight(RELATION_R4)
				pseudo_r1_freq = TaxaPair_Reln_Dict[key1]._GetFreqPseudoR1(0)
				pseudo_r2_freq = TaxaPair_Reln_Dict[key1]._GetFreqPseudoR1(1)

				if ((r1_freq + 2 * (pseudo_r1_freq - pseudo_r2_freq)) >= r4_freq) and (RELATION_R1 in allowed_reln_list):
					if (round(r1_level_val_ratio, 2) >= R1R2Reln_MAJ_THRS_very_low):	# R1R2Reln_MAJ_THRS_low):	#sourya
						return 1
					 
			if key2 in TaxaPair_Reln_Dict:
				r2_level_val_ratio = TaxaPair_Reln_Dict[key2]._GetLevelValRatio(1)
				allowed_reln_list = TaxaPair_Reln_Dict[key2]._GetAllowedRelnList()
				r2_freq = TaxaPair_Reln_Dict[key2]._GetEdgeWeight(RELATION_R2)
				r4_freq = TaxaPair_Reln_Dict[key2]._GetEdgeWeight(RELATION_R4)
				pseudo_r1_freq = TaxaPair_Reln_Dict[key2]._GetFreqPseudoR1(0)
				pseudo_r2_freq = TaxaPair_Reln_Dict[key2]._GetFreqPseudoR1(1)

				if ((r2_freq + 2 * (pseudo_r2_freq - pseudo_r1_freq)) >= r4_freq) and (RELATION_R2 in allowed_reln_list):
					if (round(r2_level_val_ratio, 2) >= R1R2Reln_MAJ_THRS_very_low):	# R1R2Reln_MAJ_THRS_low):	#sourya
						return 1

	return 0

#-------------------------------------------------------
"""
this function checks whether a cluster pair can be related by R1 or R2 relation
even if R4 relation is predominant among them
provided no conflict is induced
Here the source taxa cluster has cardinality > 1
"""
def CheckHiddenR1R2Reln(Reachability_Graph_Mat, src_taxa_clust_idx, dest_taxa_clust_idx, Output_Text_File):
	"""
	find the taxa list of respective clusters
	"""
	src_cluster_taxa_list = Cluster_Info_Dict[src_taxa_clust_idx]._GetSpeciesList()
	dest_cluster_taxa_list = Cluster_Info_Dict[dest_taxa_clust_idx]._GetSpeciesList()

	"""
	case A - src cluster size > 1
	dest cluster size = 1
	R1 relation from src cluster to dest cluster is sought
	"""
	if (len(src_cluster_taxa_list) > 1) and (len(dest_cluster_taxa_list) == 1):
		res = CheckR1Reln(src_cluster_taxa_list, dest_cluster_taxa_list)
		if (DEBUG_LEVEL >= 2):
			fp = open(Output_Text_File, 'a')
			fp.write('\n ---- src_cluster_taxa_list: ' + str(src_cluster_taxa_list) + ' dest_cluster_taxa_list: ' \
				+ str(dest_cluster_taxa_list) + '  function CheckR1Reln -- res: ' + str(res))
			fp.close()
		if (res == 1):
			if (CheckTransitiveConflict(src_taxa_clust_idx, dest_taxa_clust_idx, \
							Reachability_Graph_Mat, RELATION_R1, Output_Text_File) == 0):
				return 1, RELATION_R1
		elif (res == 2):
			Cluster_Info_Dict[src_taxa_clust_idx]._AddPossibleR1(dest_taxa_clust_idx)
			if src_taxa_clust_idx not in Candidate_Out_Edge_Cluster_List:
				Candidate_Out_Edge_Cluster_List.append(src_taxa_clust_idx)
			return 0, RELATION_R4

	"""
	case B - src cluster size = 1
	dest cluster size > 1
	R2 relation from dest cluster to src cluster is sought
	"""
	if (len(src_cluster_taxa_list) == 1) and (len(dest_cluster_taxa_list) > 1):
		res = CheckR1Reln(dest_cluster_taxa_list, src_cluster_taxa_list)
		if (DEBUG_LEVEL >= 2):
			fp = open(Output_Text_File, 'a')
			fp.write('\n ---- src_cluster_taxa_list: ' + str(src_cluster_taxa_list) + ' dest_cluster_taxa_list: ' \
				+ str(dest_cluster_taxa_list) + '  function CheckR1Reln -- res: ' + str(res))
			fp.close()
		if (res == 1):
			if (CheckTransitiveConflict(src_taxa_clust_idx, dest_taxa_clust_idx, \
							Reachability_Graph_Mat, RELATION_R2, Output_Text_File) == 0):
				return 1, RELATION_R2
		elif (res == 2):
			Cluster_Info_Dict[dest_taxa_clust_idx]._AddPossibleR1(src_taxa_clust_idx)
			if dest_taxa_clust_idx not in Candidate_Out_Edge_Cluster_List:
				Candidate_Out_Edge_Cluster_List.append(dest_taxa_clust_idx)
			return 0, RELATION_R4

	"""
	case C - src cluster size > 1
	dest cluster size > 1
	R1 / R2 relation from source to destination cluster is sought
	"""
	if ((len(src_cluster_taxa_list) > 1) and (len(dest_cluster_taxa_list) > 1)):
		res = CheckR1Reln(src_cluster_taxa_list, dest_cluster_taxa_list)
		if (res == 1):
			if (CheckTransitiveConflict(src_taxa_clust_idx, dest_taxa_clust_idx, \
							Reachability_Graph_Mat, RELATION_R1, Output_Text_File) == 0):
				return 1, RELATION_R1
		elif (res == 2):
			Cluster_Info_Dict[src_taxa_clust_idx]._AddPossibleR1(dest_taxa_clust_idx)
			if src_taxa_clust_idx not in Candidate_Out_Edge_Cluster_List:
				Candidate_Out_Edge_Cluster_List.append(src_taxa_clust_idx)
			return 0, RELATION_R4

		res = CheckR1Reln(dest_cluster_taxa_list, src_cluster_taxa_list)
		if (res == 1):
			if (CheckTransitiveConflict(src_taxa_clust_idx, dest_taxa_clust_idx, \
							Reachability_Graph_Mat, RELATION_R2, Output_Text_File) == 0):
				return 1, RELATION_R2
		elif (res == 2):
			Cluster_Info_Dict[dest_taxa_clust_idx]._AddPossibleR1(src_taxa_clust_idx)
			if dest_taxa_clust_idx not in Candidate_Out_Edge_Cluster_List:
				Candidate_Out_Edge_Cluster_List.append(dest_taxa_clust_idx)
			return 0, RELATION_R4

	#----------------------------------------------------------
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
			return 0, RELATION_R4
		
		res = CheckCandidateR1R2Reln(dest_cluster_taxa_list, src_cluster_taxa_list)
		if (res == 1):
			Cluster_Info_Dict[dest_taxa_clust_idx]._AddPossibleR1(src_taxa_clust_idx)
			if dest_taxa_clust_idx not in Candidate_Out_Edge_Cluster_List:
				Candidate_Out_Edge_Cluster_List.append(dest_taxa_clust_idx)
			return 0, RELATION_R4

	"""
	case E - src cluster size = 1 and dest cluster size > 1
	check R1 relation from src cluster to dest cluster
	"""
	if (len(src_cluster_taxa_list) == 1) and (len(dest_cluster_taxa_list) > 1):
		res = CheckCandidateR1R2Reln(src_cluster_taxa_list, dest_cluster_taxa_list)
		if (res == 1):
			Cluster_Info_Dict[src_taxa_clust_idx]._AddPossibleR1(dest_taxa_clust_idx)
			if src_taxa_clust_idx not in Candidate_Out_Edge_Cluster_List:
				Candidate_Out_Edge_Cluster_List.append(src_taxa_clust_idx)

	"""
	case F - src cluster size > 1 and dest cluster size = 1
	check R1 relation from dest cluster to src cluster
	"""
	if (len(src_cluster_taxa_list) > 1) and (len(dest_cluster_taxa_list) == 1):
		res = CheckCandidateR1R2Reln(dest_cluster_taxa_list, src_cluster_taxa_list)
		if (res == 1):
			Cluster_Info_Dict[dest_taxa_clust_idx]._AddPossibleR1(src_taxa_clust_idx)
			if dest_taxa_clust_idx not in Candidate_Out_Edge_Cluster_List:
				Candidate_Out_Edge_Cluster_List.append(dest_taxa_clust_idx)

	return 0, RELATION_R4

#-------------------------------------------------------
"""
this function checks whether R1 relation from clust1_taxa_list to clust2_taxa_list 
can be established or not, when the consensus relation is R4
"""
def CheckR1Reln(clust1_spec_list, clust2_spec_list, allowed_reln_check=True):
	# this value will be returned
	res = 0
	# explore all taxa pairs of the cluster pair
	for t1 in clust1_spec_list:
		for t2 in clust2_spec_list:
			key1 = (t1, t2)
			key2 = (t2, t1)
			if key1 in TaxaPair_Reln_Dict:
				r1_level_val_ratio = TaxaPair_Reln_Dict[key1]._GetLevelValRatio(0)
				allowed_reln_list = TaxaPair_Reln_Dict[key1]._GetAllowedRelnList()
				r1_freq = TaxaPair_Reln_Dict[key1]._GetEdgeWeight(RELATION_R1)
				r4_freq = TaxaPair_Reln_Dict[key1]._GetEdgeWeight(RELATION_R4)
				pseudo_r1_freq = TaxaPair_Reln_Dict[key1]._GetFreqPseudoR1(0)
				pseudo_r2_freq = TaxaPair_Reln_Dict[key1]._GetFreqPseudoR1(1)
				
				if (((r1_freq + 2 * (pseudo_r1_freq - pseudo_r2_freq)) >= r4_freq) \
					or (TaxaPair_Reln_Dict[key1]._CheckTargetRelnLevelConsensus(RELATION_R1))) \
						and ((RELATION_R1 in allowed_reln_list) or (allowed_reln_check == False)):
					if (round(r1_level_val_ratio, 2) >= R1R2Reln_MAJ_THRS_low):
						res = 1
					elif (round(r1_level_val_ratio, 2) >= R1R2Reln_MAJ_THRS_very_low) and (round(r1_level_val_ratio, 2) < R1R2Reln_MAJ_THRS_low):
						if (res == 0):
							res = 2
					else:
						return 0
				else:
					return 0
				
			if key2 in TaxaPair_Reln_Dict:
				r2_level_val_ratio = TaxaPair_Reln_Dict[key2]._GetLevelValRatio(1)
				allowed_reln_list = TaxaPair_Reln_Dict[key2]._GetAllowedRelnList()
				r2_freq = TaxaPair_Reln_Dict[key2]._GetEdgeWeight(RELATION_R2)
				r4_freq = TaxaPair_Reln_Dict[key2]._GetEdgeWeight(RELATION_R4)
				pseudo_r1_freq = TaxaPair_Reln_Dict[key2]._GetFreqPseudoR1(0)
				pseudo_r2_freq = TaxaPair_Reln_Dict[key2]._GetFreqPseudoR1(1)

				if (((r2_freq + 2 * (pseudo_r2_freq - pseudo_r1_freq)) >= r4_freq) \
					or (TaxaPair_Reln_Dict[key2]._CheckTargetRelnLevelConsensus(RELATION_R2))) \
						and ((RELATION_R2 in allowed_reln_list) or (allowed_reln_check == False)):
					if (round(r2_level_val_ratio, 2) >= R1R2Reln_MAJ_THRS_low):
						res = 1
					elif (round(r2_level_val_ratio, 2) >= R1R2Reln_MAJ_THRS_very_low) and (round(r2_level_val_ratio, 2) < R1R2Reln_MAJ_THRS_low):
						if (res == 0):
							res = 2
					else:
						return 0
				else:
					return 0

	return res

#-------------------------------------------------------
"""
this function checks for cluster pairs
whether they can be further analyzed for R1 / R2 relation

return: 1 if R1 reln from clust1 to clust2 is possible
				0 if no such relation is possible
				2 if immediately directed out edge from clust1 to clust2 can be established
"""
def CheckCandidateR1R2Reln(clust1_spec_list, clust2_spec_list):
	
	for t1 in clust1_spec_list:
		for t2 in clust2_spec_list:
			key1 = (t1, t2)
			key2 = (t2, t1)
			if key1 in TaxaPair_Reln_Dict:
				r1_level_val_ratio = TaxaPair_Reln_Dict[key1]._GetLevelValRatio(0)
				allowed_reln_list = TaxaPair_Reln_Dict[key1]._GetAllowedRelnList()
				#r1_freq = TaxaPair_Reln_Dict[key1]._GetEdgeWeight(RELATION_R1)
				#r4_freq = TaxaPair_Reln_Dict[key1]._GetEdgeWeight(RELATION_R4)
				#pseudo_r1_freq = TaxaPair_Reln_Dict[key1]._GetFreqPseudoR1(0)
				#pseudo_r2_freq = TaxaPair_Reln_Dict[key1]._GetFreqPseudoR1(1)
				#r4_score = TaxaPair_Reln_Dict[key1]._GetEdgeCost_ConnReln(RELATION_R4)
				
				## add - sourya
				#if (len(clust1_spec_list) == 1) and (len(clust2_spec_list) == 1):
					## in the following condition, a strict R1 relation will be established
					#if (RELATION_R1 in allowed_reln_list):
						#if ((r1_freq + pseudo_r1_freq - pseudo_r2_freq) >= r4_freq):
							#if (TaxaPair_Reln_Dict[key1]._CheckTargetRelnLevelConsensus(RELATION_R1)):
								#if (round(r1_level_val_ratio, 2) >= R1R2Reln_MAJ_THRS_very_low):
									#return 2
					### otherwise, if the support score is positive for the relation R4, then 
					### no R1 relation is possible
					##if (r4_score > 0):
						##return 0
				## end add - sourya
				
				if (RELATION_R1 in allowed_reln_list):
					if (round(r1_level_val_ratio, 2) >= R1R2Reln_MAJ_THRS_low):
						return 1
					if (TaxaPair_Reln_Dict[key1]._CheckTargetRelnLevelConsensus(RELATION_R1)) \
						and (round(r1_level_val_ratio, 2) >= R1R2Reln_MAJ_THRS_very_low):
						return 1
				
				## add - sourya
				#if (((r1_freq + 2 * (pseudo_r1_freq - pseudo_r2_freq)) >= r4_freq) \
					#or (TaxaPair_Reln_Dict[key1]._CheckTargetRelnLevelConsensus(RELATION_R1))) and (RELATION_R1 in allowed_reln_list):
					#if (round(r1_level_val_ratio, 2) >= R1R2Reln_MAJ_THRS_low):
						#res = 1
					#else:
						#return 0
				#else:
					#return 0
				## end add - sourya


			if key2 in TaxaPair_Reln_Dict:
				r2_level_val_ratio = TaxaPair_Reln_Dict[key2]._GetLevelValRatio(1)
				allowed_reln_list = TaxaPair_Reln_Dict[key2]._GetAllowedRelnList()
				#r2_freq = TaxaPair_Reln_Dict[key2]._GetEdgeWeight(RELATION_R2)
				#r4_freq = TaxaPair_Reln_Dict[key2]._GetEdgeWeight(RELATION_R4)
				#pseudo_r1_freq = TaxaPair_Reln_Dict[key2]._GetFreqPseudoR1(0)
				#pseudo_r2_freq = TaxaPair_Reln_Dict[key2]._GetFreqPseudoR1(1)
				#r4_score = TaxaPair_Reln_Dict[key2]._GetEdgeCost_ConnReln(RELATION_R4)
				
				## add - sourya
				#if (len(clust1_spec_list) == 1) and (len(clust2_spec_list) == 1):
					## in the following condition, a strict R1 relation will be established
					#if (RELATION_R2 in allowed_reln_list):
						#if ((r2_freq + pseudo_r2_freq - pseudo_r1_freq) >= r4_freq):
							#if (TaxaPair_Reln_Dict[key2]._CheckTargetRelnLevelConsensus(RELATION_R2)):
								#if (round(r2_level_val_ratio, 2) >= R1R2Reln_MAJ_THRS_very_low):
									#return 2
					### otherwise, if the support score is positive for the relation R4, then 
					### no R1 relation is possible
					##if (r4_score > 0):
						##return 0
				## end add - sourya
				
				if (RELATION_R2 in allowed_reln_list):
					if (round(r2_level_val_ratio, 2) >= R1R2Reln_MAJ_THRS_low):
						return 1

					if (TaxaPair_Reln_Dict[key2]._CheckTargetRelnLevelConsensus(RELATION_R2)) \
						and (round(r2_level_val_ratio, 2) >= R1R2Reln_MAJ_THRS_very_low):
						return 1
				# end modify - sourya

				## add - sourya
				#if (((r2_freq + 2 * (pseudo_r2_freq - pseudo_r1_freq)) >= r4_freq) \
					#or (TaxaPair_Reln_Dict[key2]._CheckTargetRelnLevelConsensus(RELATION_R2))) and (RELATION_R2 in allowed_reln_list):
					#if (round(r2_level_val_ratio, 2) >= R1R2Reln_MAJ_THRS_low):
						#res = 1
					#else:
						#return 0
				#else:
					#return 0
				## end add - sourya
	
	return 0

  