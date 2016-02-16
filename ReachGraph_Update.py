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

		## comment - sourya
		
		## for A->B connection
		## if D><A exists
		## then establish D><B
		#for x in Cluster_Info_Dict[src_taxa_clust_idx]._GetNoEdgeList():
			#if (Reachability_Graph_Mat[CURRENT_CLUST_IDX_LIST.index(x)][dest_reach_mat_idx] == 0):
				#Connect_ClusterPair(Reachability_Graph_Mat, CURRENT_CLUST_IDX_LIST.index(x), \
					#dest_reach_mat_idx, RELATION_R4, x, dest_taxa_clust_idx)
				
		## for A->B connection
		## if D><A exists
		## then for all B->E
		## establish D><E
		#for x in Cluster_Info_Dict[src_taxa_clust_idx]._GetNoEdgeList():
			#for y in Cluster_Info_Dict[dest_taxa_clust_idx]._GetOutEdgeList():
				#if (Reachability_Graph_Mat[CURRENT_CLUST_IDX_LIST.index(x)][CURRENT_CLUST_IDX_LIST.index(y)] == 0):
					#Connect_ClusterPair(Reachability_Graph_Mat, CURRENT_CLUST_IDX_LIST.index(x), \
						#CURRENT_CLUST_IDX_LIST.index(y), RELATION_R4, x, y)  
	
		## end comment - sourya

	# comment - sourya
	
	#else:
		## construct the out neighborhood of src_cluster
		## it will contain the cluster itself and all the other clusters connected via out edges from this cluster
		#src_clust_out_neighb = []
		#src_clust_out_neighb.append(src_taxa_clust_idx)
		#src_clust_out_neighb.extend(Cluster_Info_Dict[src_taxa_clust_idx]._GetOutEdgeList())
		## construct the out neighborhood of dest cluster
		## it will contain the cluster itself and all the other clusters connected via out edges from this cluster    
		#dest_clust_out_neighb = []
		#dest_clust_out_neighb.append(dest_taxa_clust_idx)
		#dest_clust_out_neighb.extend(Cluster_Info_Dict[dest_taxa_clust_idx]._GetOutEdgeList())
		
		#for x in src_clust_out_neighb:
			#for y in dest_clust_out_neighb:
				#if (Reachability_Graph_Mat[CURRENT_CLUST_IDX_LIST.index(x)][CURRENT_CLUST_IDX_LIST.index(y)] == 0):
					#Connect_ClusterPair(Reachability_Graph_Mat, CURRENT_CLUST_IDX_LIST.index(x), \
						#CURRENT_CLUST_IDX_LIST.index(y), RELATION_R4, x, y)  
	
	# end comment - sourya
	
#-----------------------------------------------------
""" 
this function merges two clusters 
basically the cluster src_clust_idx will be merged to the dest_clust_idx
so all the entries concerning src_clust_idx will now point to the dest_clust_idx
"""
def Merge_Clusters(Reachability_Graph_Mat, dest_taxa_label, src_taxa_label, dest_clust_idx, src_clust_idx, \
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
def AdjustReachGraph(Reachability_Graph_Mat, taxaA_label, taxaB_label, reln_type):    

	nodeA_clust_idx = Taxa_Info_Dict[taxaA_label]._Get_Taxa_Part_Clust_Idx()
	nodeB_clust_idx = Taxa_Info_Dict[taxaB_label]._Get_Taxa_Part_Clust_Idx()

	# perform cluster merging operation (content shifting plus transitive closure)
	nodeA_reach_mat_idx = CURRENT_CLUST_IDX_LIST.index(nodeA_clust_idx)
	nodeB_reach_mat_idx = CURRENT_CLUST_IDX_LIST.index(nodeB_clust_idx)

	# update tge reachability matrix information
	if (reln_type == RELATION_R3):      
		#-----------------------------
		# keep the minimum cluster index intact
		# merge two clusters
		if (nodeA_clust_idx > nodeB_clust_idx):
			Merge_Clusters(Reachability_Graph_Mat, taxaB_label, taxaA_label, nodeB_clust_idx, nodeA_clust_idx, nodeB_reach_mat_idx, nodeA_reach_mat_idx)
			# delete the index of nodeA_clust_idx from the reachability matrix
			Reachability_Graph_Mat = numpy.delete(Reachability_Graph_Mat, (nodeA_reach_mat_idx), axis=0)	# delete the row
			Reachability_Graph_Mat = numpy.delete(Reachability_Graph_Mat, (nodeA_reach_mat_idx), axis=1)	# delete the column
			# delete the entry of nodeA_clust_idx from the CURRENT_CLUST_IDX_LIST
			CURRENT_CLUST_IDX_LIST.remove(nodeA_clust_idx)
			# also remove the cluster key from the dictionary
			Cluster_Info_Dict.pop(nodeA_clust_idx, None)          
		else:
			Merge_Clusters(Reachability_Graph_Mat, taxaA_label, taxaB_label, nodeA_clust_idx, nodeB_clust_idx, nodeA_reach_mat_idx, nodeB_reach_mat_idx)    
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

#-------------------------------------------------------------
# add - sourya

"""
this function processes all clusters having candidate out edge information
"""
def Process_Candidate_Out_Edge_Cluster_List(Reachability_Graph_Mat, DIST_MAT_TYPE):
	
	for cl in Candidate_Out_Edge_Cluster_List:
		"""
		check if the cluster cl has any directed out edge already
		otherwise, place all its candidate R1 clusters as a child to its parent node
		"""
		if (Cluster_Info_Dict[cl]._Get_Outdegree() == 0):
			"""
			here, assign out edge from the parent cluster of cl to the candidate R1 clusters
			"""
			for parent_cl in Cluster_Info_Dict[cl]._GetInEdgeList():
				parent_cl_reach_mat_idx = CURRENT_CLUST_IDX_LIST.index(parent_cl)
				for x in Cluster_Info_Dict[cl]._GetPossibleR1List():
					x_reach_mat_idx = CURRENT_CLUST_IDX_LIST.index(x)
					# connect the pair of clusters, along with updating the reachability matrix
					Connect_ClusterPair(Reachability_Graph_Mat, parent_cl_reach_mat_idx, x_reach_mat_idx, RELATION_R1, parent_cl, x)
					# now perform the transitive closure on the derived reachability matrix  
					TransClosUpd(Reachability_Graph_Mat, parent_cl_reach_mat_idx, x_reach_mat_idx, RELATION_R1)
		else:
			"""
			explore children of the cluster cl
			and compare its XL with the new cluster
			"""
			cl_reach_mat_idx = CURRENT_CLUST_IDX_LIST.index(cl)
			
			for x in Cluster_Info_Dict[cl]._GetPossibleR1List():
				x_reach_mat_idx = CURRENT_CLUST_IDX_LIST.index(x)
				"""
				excess gene count between the cluster x and the cluster cl
				"""
				xl_cl_x = FindAvgXL(Cluster_Info_Dict[cl]._GetSpeciesList(), Cluster_Info_Dict[x]._GetSpeciesList(), DIST_MAT_TYPE, True)
				"""
				this is the average of excess gene count between x and every child of cl
				"""
				xl_x_childcl = 0
				"""
				this is the average of excess gene count between cl and every child of cl
				"""
				xl_cl_childcl = 0
				
				for child_cl in Cluster_Info_Dict[cl]._GetOutEdgeList():
					xl_x_childcl = xl_x_childcl + FindAvgXL(Cluster_Info_Dict[child_cl]._GetSpeciesList(), \
						Cluster_Info_Dict[x]._GetSpeciesList(), DIST_MAT_TYPE, True)
					xl_cl_childcl = xl_cl_childcl + FindAvgXL(Cluster_Info_Dict[child_cl]._GetSpeciesList(), \
						Cluster_Info_Dict[cl]._GetSpeciesList(), DIST_MAT_TYPE, True)
				
				xl_x_childcl = (xl_x_childcl * 1.0) / len(Cluster_Info_Dict[cl]._GetOutEdgeList())
				xl_cl_childcl = (xl_cl_childcl * 1.0) / len(Cluster_Info_Dict[cl]._GetOutEdgeList())
				
				if (xl_x_childcl <= xl_cl_x) and (xl_x_childcl <= xl_cl_childcl):
					"""
					x can be placed as the child of cl
					"""
					# connect the pair of clusters, along with updating the reachability matrix
					Connect_ClusterPair(Reachability_Graph_Mat, cl_reach_mat_idx, x_reach_mat_idx, RELATION_R1, cl, x)
					# now perform the transitive closure on the derived reachability matrix  
					TransClosUpd(Reachability_Graph_Mat, cl_reach_mat_idx, x_reach_mat_idx, RELATION_R1)
				else:
					"""
					here, assign out edge from the parent cluster of cl to the candidate R1 clusters
					"""
					for parent_cl in Cluster_Info_Dict[cl]._GetInEdgeList():
						parent_cl_reach_mat_idx = CURRENT_CLUST_IDX_LIST.index(parent_cl)
						for x in Cluster_Info_Dict[cl]._GetPossibleR1List():
							x_reach_mat_idx = CURRENT_CLUST_IDX_LIST.index(x)
							# connect the pair of clusters, along with updating the reachability matrix
							Connect_ClusterPair(Reachability_Graph_Mat, parent_cl_reach_mat_idx, x_reach_mat_idx, RELATION_R1, parent_cl, x)
							# now perform the transitive closure on the derived reachability matrix  
							TransClosUpd(Reachability_Graph_Mat, parent_cl_reach_mat_idx, x_reach_mat_idx, RELATION_R1)
					
					
				