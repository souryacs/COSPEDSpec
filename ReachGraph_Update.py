#!/usr/bin/env python

import Header
from Header import *
import Cluster_Manage
from Cluster_Manage import *

#-----------------------------------------------------
""" 
this function adds an edge between a pair of clusters (of taxa) 
it also updates the entries of reachability matrix 
"""
def Connect_ClusterPair(Reachability_Graph_Mat, nodeA_reach_mat_idx, nodeB_reach_mat_idx, reln_type, nodeA_clust_idx, nodeB_clust_idx):
	if (reln_type == RELATION_R1):
		# adjust the clusters
		Cluster_Info_Dict[nodeA_clust_idx]._AddOutEdge(nodeB_clust_idx)
		Cluster_Info_Dict[nodeB_clust_idx]._AddInEdge(nodeA_clust_idx)
		# update the reachability matrix
		Reachability_Graph_Mat[nodeA_reach_mat_idx][nodeB_reach_mat_idx] = 1
	elif (reln_type == RELATION_R4):
		# adjust the clusters
		Cluster_Info_Dict[nodeA_clust_idx]._AddNoEdge(nodeB_clust_idx)
		Cluster_Info_Dict[nodeB_clust_idx]._AddNoEdge(nodeA_clust_idx)    
		# update the reachability matrix
		Reachability_Graph_Mat[nodeA_reach_mat_idx][nodeB_reach_mat_idx] = 2
		Reachability_Graph_Mat[nodeB_reach_mat_idx][nodeA_reach_mat_idx] = 2
        
#-----------------------------------------------------
""" this function updates the transitive closure of the cluster of nodes
on inclusion ogf a new edge between a pair of clusters """
def TransClosUpd(Reachability_Graph_Mat, nodeA_reach_mat_idx, nodeB_reach_mat_idx, reln_type):
	if (reln_type == RELATION_R1) or (reln_type == RELATION_R4):
		src_reach_mat_idx = nodeA_reach_mat_idx
		dest_reach_mat_idx = nodeB_reach_mat_idx
	elif (reln_type == RELATION_R2):
		src_reach_mat_idx = nodeB_reach_mat_idx
		dest_reach_mat_idx = nodeA_reach_mat_idx
	else:
		return
		
	src_taxa_clust_idx = CURRENT_CLUST_IDX_LIST[src_reach_mat_idx]
	dest_taxa_clust_idx = CURRENT_CLUST_IDX_LIST[dest_reach_mat_idx]
		
	if (reln_type == RELATION_R1) or (reln_type == RELATION_R2):
		# for A->B connection
		# if D->A exists
		# then establish D->B
		for x in Cluster_Info_Dict[src_taxa_clust_idx]._GetInEdgeList():
			if (Reachability_Graph_Mat[CURRENT_CLUST_IDX_LIST.index(x)][dest_reach_mat_idx] == 0):
				Connect_ClusterPair(Reachability_Graph_Mat, CURRENT_CLUST_IDX_LIST.index(x), \
					dest_reach_mat_idx, RELATION_R1, x, dest_taxa_clust_idx)
		
		# for A->B connection
		# if B->E exists
		# then establish A->E
		for x in Cluster_Info_Dict[dest_taxa_clust_idx]._GetOutEdgeList():
			if (Reachability_Graph_Mat[src_reach_mat_idx][CURRENT_CLUST_IDX_LIST.index(x)] == 0):
				Connect_ClusterPair(Reachability_Graph_Mat, src_reach_mat_idx, \
					CURRENT_CLUST_IDX_LIST.index(x), RELATION_R1, src_taxa_clust_idx, x)

		# for A->B connection
		# if D->A and B->E exists
		# then establish D->E  
		for x in Cluster_Info_Dict[src_taxa_clust_idx]._GetInEdgeList():
			for y in Cluster_Info_Dict[dest_taxa_clust_idx]._GetOutEdgeList():
				if (Reachability_Graph_Mat[CURRENT_CLUST_IDX_LIST.index(x)][CURRENT_CLUST_IDX_LIST.index(y)] == 0):  
					Connect_ClusterPair(Reachability_Graph_Mat, \
						CURRENT_CLUST_IDX_LIST.index(x), CURRENT_CLUST_IDX_LIST.index(y), RELATION_R1, x, y)  

#-----------------------------------------------------
""" 
this function merges two clusters 
basically the cluster src_clust_idx will be merged to the dest_clust_idx
so all the entries concerning src_clust_idx will now point to the dest_clust_idx
"""
def Merge_Clusters(Reachability_Graph_Mat, dest_clust_idx, src_clust_idx, \
		    dest_clust_reach_mat_idx, src_clust_reach_mat_idx):
	
	""" 
	first update the reachability matrix entries 
	originally all the Reachability_Graph_Mat entries concerning src_clust_reach_mat_idx (for the cluster src_clust_idx)
	will now point to the dest_clust_reach_mat_idx (for the cluster dest_clust_idx) as well
	also update the dest cluster out edge, in edge, and no edge lists (corresponding to the relations R1, R2, and R4)
	"""
	
	"""
	Important - Note -
	In copying the out / in / no edge information
	we can overwrite the no edge with a definite out / in edge
	"""

	for x in Cluster_Info_Dict[src_clust_idx]._GetOutEdgeList():
		Cluster_Info_Dict[x]._RemoveInEdge(src_clust_idx)
		#Cluster_Info_Dict[src_clust_idx]._RemoveOutEdge(x)
		if x in Cluster_Info_Dict[dest_clust_idx]._GetNoEdgeList():
			Cluster_Info_Dict[x]._RemoveNoEdge(dest_clust_idx)
			Cluster_Info_Dict[dest_clust_idx]._RemoveNoEdge(x)
		if (Reachability_Graph_Mat[dest_clust_reach_mat_idx][CURRENT_CLUST_IDX_LIST.index(x)] != 1):
			Connect_ClusterPair(Reachability_Graph_Mat, dest_clust_reach_mat_idx, CURRENT_CLUST_IDX_LIST.index(x), RELATION_R1, dest_clust_idx, x)  

	for x in Cluster_Info_Dict[src_clust_idx]._GetInEdgeList():
		Cluster_Info_Dict[x]._RemoveOutEdge(src_clust_idx)
		#Cluster_Info_Dict[src_clust_idx]._RemoveInEdge(x)
		if x in Cluster_Info_Dict[dest_clust_idx]._GetNoEdgeList():
			Cluster_Info_Dict[x]._RemoveNoEdge(dest_clust_idx)
			Cluster_Info_Dict[dest_clust_idx]._RemoveNoEdge(x)
		if (Reachability_Graph_Mat[CURRENT_CLUST_IDX_LIST.index(x)][dest_clust_reach_mat_idx] != 1):
			Connect_ClusterPair(Reachability_Graph_Mat, CURRENT_CLUST_IDX_LIST.index(x), dest_clust_reach_mat_idx, RELATION_R1, x, dest_clust_idx)  

	for x in Cluster_Info_Dict[src_clust_idx]._GetNoEdgeList():
		Cluster_Info_Dict[x]._RemoveNoEdge(src_clust_idx)
		#Cluster_Info_Dict[src_clust_idx]._RemoveNoEdge(x)
		if (Reachability_Graph_Mat[CURRENT_CLUST_IDX_LIST.index(x)][dest_clust_reach_mat_idx] == 0):
			Connect_ClusterPair(Reachability_Graph_Mat, CURRENT_CLUST_IDX_LIST.index(x), dest_clust_reach_mat_idx, RELATION_R4, x, dest_clust_idx)      
	
	# then adjust the taxa in the other cluster
	for tax in Cluster_Info_Dict[src_clust_idx]._GetSpeciesList():
		Append_Cluster_Taxa_Label(dest_clust_idx, tax)
	
#-----------------------------------------------------
""" this function updates the reachability graph 
on the basis of input edge type between input 2 taxa """
def AdjustReachGraph(Reachability_Graph_Mat, nodeA_clust_idx, nodeB_clust_idx, reln_type):    

	#nodeA_clust_idx = Taxa_Info_Dict[taxaA_label]._Get_Taxa_Part_Clust_Idx()
	#nodeB_clust_idx = Taxa_Info_Dict[taxaB_label]._Get_Taxa_Part_Clust_Idx()

	# perform cluster merging operation (content shifting plus transitive closure)
	nodeA_reach_mat_idx = CURRENT_CLUST_IDX_LIST.index(nodeA_clust_idx)
	nodeB_reach_mat_idx = CURRENT_CLUST_IDX_LIST.index(nodeB_clust_idx)

	# update tge reachability matrix information
	if (reln_type == RELATION_R3):      
		#-----------------------------
		# keep the minimum cluster index intact
		# merge two clusters
		if (nodeA_clust_idx > nodeB_clust_idx):
			Merge_Clusters(Reachability_Graph_Mat, nodeB_clust_idx, nodeA_clust_idx, nodeB_reach_mat_idx, nodeA_reach_mat_idx)
			# delete the index of nodeA_clust_idx from the reachability matrix
			Reachability_Graph_Mat = numpy.delete(Reachability_Graph_Mat, (nodeA_reach_mat_idx), axis=0)	# delete the row
			Reachability_Graph_Mat = numpy.delete(Reachability_Graph_Mat, (nodeA_reach_mat_idx), axis=1)	# delete the column
			# delete the entry of nodeA_clust_idx from the CURRENT_CLUST_IDX_LIST
			CURRENT_CLUST_IDX_LIST.remove(nodeA_clust_idx)
			# also remove the cluster key from the dictionary
			Cluster_Info_Dict.pop(nodeA_clust_idx, None)          
		else:
			Merge_Clusters(Reachability_Graph_Mat, nodeA_clust_idx, nodeB_clust_idx, nodeA_reach_mat_idx, nodeB_reach_mat_idx)    
			# delete the index of nodeB_clust_idx from the reachability matrix
			Reachability_Graph_Mat = numpy.delete(Reachability_Graph_Mat, (nodeB_reach_mat_idx), axis=0)	# delete the row
			Reachability_Graph_Mat = numpy.delete(Reachability_Graph_Mat, (nodeB_reach_mat_idx), axis=1)	# delete the column
			# delete the entry of nodeB_clust_idx from the CURRENT_CLUST_IDX_LIST
			CURRENT_CLUST_IDX_LIST.remove(nodeB_clust_idx)
			# also remove the cluster key from the dictionary
			Cluster_Info_Dict.pop(nodeB_clust_idx, None)    
		#-----------------------------
	elif (reln_type == RELATION_R1):
		# connect the pair of clusters, along with updating the reachability matrix
		Connect_ClusterPair(Reachability_Graph_Mat, nodeA_reach_mat_idx, nodeB_reach_mat_idx, \
			RELATION_R1, nodeA_clust_idx, nodeB_clust_idx)
		# now perform the transitive closure on the derived reachability matrix  
		TransClosUpd(Reachability_Graph_Mat, nodeA_reach_mat_idx, nodeB_reach_mat_idx, RELATION_R1)
	elif (reln_type == RELATION_R2):
		# connect the pair of clusters, along with updating the reachability matrix
		Connect_ClusterPair(Reachability_Graph_Mat, nodeB_reach_mat_idx, nodeA_reach_mat_idx, \
			RELATION_R1, nodeB_clust_idx, nodeA_clust_idx)
		# now perform the transitive closure on the derived reachability matrix  
		TransClosUpd(Reachability_Graph_Mat, nodeA_reach_mat_idx, nodeB_reach_mat_idx, RELATION_R2)    
	else:	#reln_type == RELATION_R4:
		Connect_ClusterPair(Reachability_Graph_Mat, nodeA_reach_mat_idx, nodeB_reach_mat_idx, \
			RELATION_R4, nodeA_clust_idx, nodeB_clust_idx)
		# now perform the transitive closure on the derived reachability matrix  
		TransClosUpd(Reachability_Graph_Mat, nodeA_reach_mat_idx, nodeB_reach_mat_idx, RELATION_R4)    
			
	return Reachability_Graph_Mat

