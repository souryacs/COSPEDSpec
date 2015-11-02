#!/usr/bin/env python

import dendropy
from dendropy import Tree, Taxon, TaxonSet, Node
import numpy
import time
import os
from cStringIO import StringIO
from optparse import OptionParser
import math
import sys

# we define custom relation types
RELATION_R3 = 0	# relation r3
RELATION_R1 = 1	# relation r1
RELATION_R2 = 2	# relation r2
RELATION_R4 = 3	# relation r4
UNDEFINED_EDGE = 4

#----------------------------------------------------------
""" 
this is a dictionary for storing information about individual taxa clusters
each cluster is basically a collection of taxa related via relation r3
"""
Cluster_Info_Dict = dict()

""" 
the dictionary defines one particular taxa and its associated information
"""
Taxa_Info_Dict = dict()

""" 
this dictionary defines the taxa pair (couplet) relations and associated operations
each entry of this dictionary is indexed by a pair of taxon labels 
"""
TaxaPair_Reln_Dict = dict()

"""
queue storing relations of conflicting couplets
"""
Cost_List_Taxa_Pair_Multi_Reln = []

""" 
queue storing relations of non-conflicting couplets
provided that we use this queue
"""
Cost_List_Taxa_Pair_Single_Reln = [] 

""" 
this list contains the complete set of taxa present in the input trees 
"""
COMPLETE_INPUT_TAXA_LIST = []

""" 
this list contains the current set of active taxa cluster (indices)
"""
CURRENT_CLUST_IDX_LIST = []

# this is the debug level
# set for printing the necessary information
DEBUG_LEVEL = 0

# this text file stores all the printing output
Output_Text_File = 'complete_output_description.txt'

FREQ_COUNT_PERCENT_THR_R4 = 0.55
FREQ_COUNT_PERCENT_THR_R1R2 = 0.45
PERCENT_R1R2_LEV_COUNT = 0.8

##-----------------------------------------------------
""" 
this class defines a taxon
"""
class Single_Taxa(object):
	def __init__(self):
		""" 
		this variable signifies the cluster index that the taxa belongs 
		if the value is -1 (< 0) then it is still not inserted in a cluster 
		otherwise it is part of a valid cluster 
		"""
		self.clust_idx_part = -1

	def _Get_Taxa_Part_Clust_Idx(self):
		return self.clust_idx_part
		
	def _Set_Clust_Idx_taxa_Part(self, inp_clust_idx):
		self.clust_idx_part = inp_clust_idx
			
	# this function is called after formation of consensus tree
	def _PrintFinalTaxaInfo(self, key, Output_Text_File):
		fp = open(Output_Text_File, 'a')    
		fp.write('\n taxa key: ' + str(key))
		fp.write('\n taxa is part of the cluster ID: ' + str(self.clust_idx_part))
		fp.close()
        
##-----------------------------------------------------
""" 
this class defines a couplet, according to the information obtained from input trees
key of this class --- taxa1, taxa2  
the consensus relation, frequencies, priority, and support scores 
in the class, the edge type signifies the relation between a pair of taxa
"""
class Reln_TaxaPair(object):
	def __init__(self):    
		
		""" 
		frequencies of individual relations 
		there are 4 types of edges (relationship) between a pair of taxa 
		"""
		self.freq_count = [0] * 4    
		
		""" 
		a connection priority value is defined as the 
		no of occurrences of this particular relation between this pair of taxa 
		minus the sum of no of occurrences of other relation types between this couplet
		"""
		self.priority_reln = [0] * 4    
		
		""" 
		this is the support score for different types of relations between a couplet
		"""
		self.support_score = [0] * 4
		
		""" 
		this variable stores the no of trees supporting the taxa pair 
		"""
		self.supporting_trees = 0
		
		""" 
		For this couplet, it stores the extra gene count with respect to all the gene trees
		"""
		self.XL_sum_gene_trees = 0
		
		"""
		for a couplet xy, and for a relation r3, following array has 3 elements:
		1) count when level of x < level of y (relation r1 similar)
		2) count when level of x > level of y (relation r2 similar)
		3) count when level of x = level of y (relation r3 similar)
		"""
		#self.R4_Reln_Level_Diff_Info_Count = [0] * 3
		#self.R4_Reln_Level_Diff_Val_Count = [0] * 3
		
		#self.ALL_Reln_Level_Diff_Info_Count = [0] * 3
		#self.ALL_Reln_Level_Diff_Val_Count = [0] * 3
		
	#def _GetAllRelnLevelDiffCount(self):
		#return self.ALL_Reln_Level_Diff_Info_Count
		
	#def _GetAllRelnLevelDiffVal(self):
		#return self.ALL_Reln_Level_Diff_Val_Count
		
	#def _IncrAllRelnLevelDiffInfoCount(self, idx, val):
		#self.ALL_Reln_Level_Diff_Info_Count[idx] = self.ALL_Reln_Level_Diff_Info_Count[idx] + 1
		#self.ALL_Reln_Level_Diff_Val_Count[idx] = self.ALL_Reln_Level_Diff_Val_Count[idx] + val
		
	#def _IncrLevelDiffInfoCount(self, idx, val):
		#self.R4_Reln_Level_Diff_Info_Count[idx] = self.R4_Reln_Level_Diff_Info_Count[idx] + 1
		#self.R4_Reln_Level_Diff_Val_Count[idx] = self.R4_Reln_Level_Diff_Val_Count[idx] + val
		
	def _AddXLVal(self, XL_val):
		self.XL_sum_gene_trees = self.XL_sum_gene_trees + XL_val
		
	def _GetXLSumGeneTrees(self):
		return self.XL_sum_gene_trees

	def _GetNormalizedXLSumGeneTrees(self):
		return (self.XL_sum_gene_trees * 1.0) / self.supporting_trees
		
	def _AddSupportingTree(self):
		self.supporting_trees = self.supporting_trees + 1
				
	def _GetEdgeWeight(self, reln_type):
		return self.freq_count[reln_type]      
		
	def _GetEdgeCost_ConnReln(self, reln_type):
		return self.support_score[reln_type]
		
	def _IncrEdgeCost_ConnReln(self, reln_type, incr_cost):
		self.support_score[reln_type] = self.support_score[reln_type] + incr_cost

	# this function adds one frequency count (with a given input relation type)
	def _AddEdgeCount(self, reln_type):
		self.freq_count[reln_type] = self.freq_count[reln_type] + 1
		
	# this function prints the relationship information
	def _PrintRelnInfo(self, key, Output_Text_File):
		fp = open(Output_Text_File, 'a')    
		fp.write('\n taxa pair key: ' + str(key))
		fp.write('\n relations [type/count/priority_reln/score]: ')
		for i in range(4):
			fp.write('\n [' + str(i) + '/' + str(self.freq_count[i]) + '/' + str(self.priority_reln[i]) + '/' + str(self.support_score[i]) + ']')
		fp.write('\n Sum of extra lineage : ' + str(self.XL_sum_gene_trees))
		fp.write('\n No of supporting trees : ' + str(self.supporting_trees))
		fp.write('\n Normalized XL sum : ' + str(self._GetNormalizedXLSumGeneTrees()))
		#fp.write('\n R4 relation based Level diff info count (r1/r2/r3): ' + str(self.R4_Reln_Level_Diff_Info_Count))
		#fp.write('\n R4 relation based Level diff Val count (r1/r2/r3): ' + str(self.R4_Reln_Level_Diff_Val_Count))
		#fp.write('\n ALL relation based Level diff info count (r1/r2/r3): ' + str(self.ALL_Reln_Level_Diff_Info_Count))
		#fp.write('\n ALL relation based Level diff Val count (r1/r2/r3): ' + str(self.ALL_Reln_Level_Diff_Val_Count))
		fp.close()
					
	# this function computes the support score metric value associated with individual pair of taxa 
	def _SetCostMetric(self):      
		for reln_type in range(4):
			# assign the score metric for this edge type
			self.support_score[reln_type] = self.freq_count[reln_type] * self.priority_reln[reln_type]
			#self.support_score[reln_type] = self.freq_count[reln_type]	# sourya - debug
			
	# this function returns the connection priority value for input relation 
	def _GetConnPrVal(self, reln_type):
		return self.priority_reln[reln_type]
	
	#"""
	#this function checks whether R4 relation is the predominant / consensus
	#in such a case, it applies level count analysis to redistribute the frequency measures
	#"""
	#def _AdjustFreq(self, key, Output_Text_File):
		##fp = open(Output_Text_File, 'a')    
		##fp.write('\n taxa pair key: ' + str(key))
		#max_freq = max(self.freq_count)
		#if (self.freq_count[RELATION_R4] == max_freq):
			## here R4 is the consensus relation
			#if (DEBUG_LEVEL >= 2):
				#fp = open(Output_Text_File, 'a')
				#fp.write('\n Couplet key: ' + str(key))
				#fp.write('\n consensus relation: ' + str(RELATION_R4))
				#fp.close()
			##------------------------------------------------
			#if (self.priority_reln[RELATION_R4] <= 0):
				#if (DEBUG_LEVEL >= 2):
					#fp = open(Output_Text_File, 'a')
					#fp.write('\n this relation is not a majority consensus relation - here R4 priority is <= 0 ')
					#fp.close()
				#""" 
				#current relation has priority less than or equal to zero
				#so this relation is although consensus, it is not majority consensus
				#now we have to check the R4_Reln_Level_Diff_Info_Count values
				#"""
				#if (self.R4_Reln_Level_Diff_Info_Count[0] > (self.R4_Reln_Level_Diff_Info_Count[1] + self.R4_Reln_Level_Diff_Info_Count[2])) \
					#or (self.R4_Reln_Level_Diff_Info_Count[1] > (self.R4_Reln_Level_Diff_Info_Count[0] + self.R4_Reln_Level_Diff_Info_Count[2])):
					#if (DEBUG_LEVEL >= 2):
						#fp = open(Output_Text_File, 'a')
						#fp.write('\n ---- Current couplet frequency distribution before modification --- ')
						#fp.close()
						#self._PrintRelnInfo(key, Output_Text_File)
					#"""
					#level count corresponding to the relation r1 / r2 is predominant
					#so move the r4 relation instance to both r1 and r2 relation instances
					#"""
					#self.freq_count[RELATION_R1] = self.freq_count[RELATION_R1] + self.R4_Reln_Level_Diff_Info_Count[0]
					#self.freq_count[RELATION_R2] = self.freq_count[RELATION_R2] + self.R4_Reln_Level_Diff_Info_Count[1]
					#self.freq_count[RELATION_R4] = self.freq_count[RELATION_R4] - self.R4_Reln_Level_Diff_Info_Count[0] - self.R4_Reln_Level_Diff_Info_Count[1]
					## recompute the priority
					#self._SetConnPrVal()
					#if (DEBUG_LEVEL >= 2):
						#fp = open(Output_Text_File, 'a')
						#fp.write('\n ---- Current couplet frequency distribution after modification --- ')
						#fp.close()
						#self._PrintRelnInfo(key, Output_Text_File)
			#else:
				#total_freq = sum(self.freq_count)
				#r1_freq = self.freq_count[RELATION_R1]
				#r2_freq = self.freq_count[RELATION_R2]
				#r3_freq = self.freq_count[RELATION_R3]
				#r4_freq = self.freq_count[RELATION_R4]
				
				#r1_lev_count = self.R4_Reln_Level_Diff_Info_Count[0]
				#r2_lev_count = self.R4_Reln_Level_Diff_Info_Count[1]
				#r4_lev_count = self.R4_Reln_Level_Diff_Info_Count[2]
				
				#sum_level_diff_val = sum(self.R4_Reln_Level_Diff_Val_Count)
				#r1_level_diff_val = self.R4_Reln_Level_Diff_Val_Count[0]
				#r2_level_diff_val = self.R4_Reln_Level_Diff_Val_Count[1]
				
				## R4 is a consensus relation but its priority is not that high
				#if (r4_freq <= (FREQ_COUNT_PERCENT_THR_R4 * total_freq)):
					#if (DEBUG_LEVEL >= 2):
						#fp = open(Output_Text_File, 'a')
						#fp.write('\n ---- R4 is a majority consensus relation but with low priority --- ')
						#fp.close()
					## case 1 - count of R1 is quite close
					#if (r1_freq >= (FREQ_COUNT_PERCENT_THR_R1R2 * total_freq)):	# and (r2_freq == 0) and (r3_freq == 0):
						## also the couplet exhibits R1 relation like clade for most of the gene trees
						#if (r1_lev_count > (r2_lev_count + r4_lev_count)) and (r1_level_diff_val >= PERCENT_R1R2_LEV_COUNT * sum_level_diff_val):
						##if (r1_level_diff_val >= PERCENT_R1R2_LEV_COUNT * sum_level_diff_val):
							## further condition
							#if ((r1_freq + r1_lev_count - r2_lev_count - r4_lev_count) > r4_freq):
								#if (DEBUG_LEVEL >= 2):
									#fp = open(Output_Text_File, 'a')
									#fp.write('\n ---- Current couplet frequency distribution before modification --- ')
									#fp.close()
									#self._PrintRelnInfo(key, Output_Text_File)
								
								#self.freq_count[RELATION_R1] = r1_freq + r1_lev_count
								#self.freq_count[RELATION_R2] = r2_freq + r2_lev_count
								#self.freq_count[RELATION_R4] = r4_freq - r1_lev_count - r2_lev_count
								## recompute the priority
								#self._SetConnPrVal()
								
								#if (DEBUG_LEVEL >= 2):
									#fp = open(Output_Text_File, 'a')
									#fp.write('\n ---- Current couplet frequency distribution after modification --- ')
									#fp.close()
									#self._PrintRelnInfo(key, Output_Text_File)

					## case 2 - count of R2 is quite close
					#if (r2_freq >= (FREQ_COUNT_PERCENT_THR_R1R2 * total_freq)):	# and (r1_freq == 0) and (r3_freq == 0):
						## also the couplet exhibits R1 relation like clade for most of the gene trees
						#if (r2_lev_count > (r1_lev_count + r4_lev_count)) and (r2_level_diff_val >= PERCENT_R1R2_LEV_COUNT * sum_level_diff_val):
						##if (r2_level_diff_val >= PERCENT_R1R2_LEV_COUNT * sum_level_diff_val):
							## further condition
							#if ((r2_freq + r2_lev_count - r1_lev_count - r4_lev_count) > r4_freq):
								#if (DEBUG_LEVEL >= 2):
									#fp = open(Output_Text_File, 'a')
									#fp.write('\n ---- Current couplet frequency distribution before modification --- ')
									#fp.close()
									#self._PrintRelnInfo(key, Output_Text_File)
								
								#self.freq_count[RELATION_R1] = r1_freq + self.R4_Reln_Level_Diff_Info_Count[0]
								#self.freq_count[RELATION_R2] = r2_freq + r2_lev_count
								#self.freq_count[RELATION_R4] = r4_freq - r1_lev_count - r2_lev_count
								## recompute the priority
								#self._SetConnPrVal()
								
								#if (DEBUG_LEVEL >= 2):
									#fp = open(Output_Text_File, 'a')
									#fp.write('\n ---- Current couplet frequency distribution after modification --- ')
									#fp.close()
									#self._PrintRelnInfo(key, Output_Text_File)
								
	""" 
	this function calculates connection priority value for each of the relation types, 
	"""
	def _SetConnPrVal(self):
		# this is the sum of frequencies for all the relation types
		listsum = sum(self.freq_count)
		# now determine the connection priority of a particular relation type with respect to other relations     
		for reln_type in range(4):
			# here we use the difference of current relation type frequency with the frequencies of all other relations
			self.priority_reln[reln_type] = 2 * self.freq_count[reln_type] - listsum
			
	def _CheckNonConflictingCouplet(self):
		# this is the sum of frequencies for all the relation types
		listsum = sum(self.freq_count)
		""" 
		this code section is used when there exists an unique relation (non-conflicting couplet)
		and we try to detect it
		"""
		outlist = [0, RELATION_R4]
		for reln_type in range(4):
			if (self.freq_count[reln_type] == listsum) and (listsum > 0):
				outlist = [1, reln_type]
				break
			elif (self.freq_count[reln_type] > 0) and (self.freq_count[reln_type] < listsum):
				break
		return outlist

##-----------------------------------------------------
""" this class is representative of a cluster of taxa that are related via equality relationship 
according to the rule of equivalence partition """
class Cluster_node(object):
	def __init__(self, inp_taxa=None):
		# this list contains the species list of the current cluster
		self.Species_List = [] 
		# can be 0 or 1 - denote whether the cluster node has been explored during newick string construction
		self.explored = 0    
		# this variable stores the out edge list from this cluster
		# each list element is the other cluster index
		self.out_edge_list = []
		# this variable stores the in edge list from this cluster
		# each list element is the other cluster index 
		self.in_edge_list = []
		# this variable stores the NO edge list from this cluster
		# each list element is the other cluster index 
		self.no_edge_list = []
		# during initialization, append one tuple to this cluster
		if inp_taxa is not None:
			self._Append_taxa(inp_taxa)    

	# these functions keep track whether the cluster node is used during newick string formation for supertree construction
	# each of the clusters (containing a set of taxa) should be visited exactly once for supertree generation
	def _SetExploredStatus(self):
		self.explored = 1

	def _ResetExploredStatus(self):
		self.explored = 0
		
	def _GetExploredStatus(self):
		return self.explored
					
	# returns the constituent species list of this cluster
	def _GetSpeciesList(self):
		return self.Species_List
				
	# append one species information in this cluster
	def _Append_taxa(self, inp_taxa):
		if inp_taxa not in self.Species_List:
			self.Species_List.append(inp_taxa)

	# it returns the final cluster node connectivity (tree shape) -- in edges
	def _Get_Indegree(self):
		return len(self.in_edge_list)

	# it returns the final cluster node connectivity (tree shape) -- out edges
	def _Get_Outdegree(self):
		return len(self.out_edge_list)
			
	# it returns the out edge -- of the cluster node to the other nodes (clique formation)
	def _GetOutEdgeList(self):
		return self.out_edge_list
		
	# it returns the in edge -- of the cluster node to the other nodes (clique formation)
	def _GetInEdgeList(self):
		return self.in_edge_list    

	# it returns the in edge -- of the cluster node to the other nodes (clique formation)
	def _GetNoEdgeList(self):
		return self.no_edge_list    
		
	# it adds one out edge information to both the original connectivity (clique) and the final connectivity (tree shape)
	def _AddOutEdge(self, dest_clust_idx):
		if dest_clust_idx not in self.out_edge_list:
			self.out_edge_list.append(dest_clust_idx)
		
	# it adds one in edge information to both the original connectivity (clique) and the final connectivity (tree shape)
	def _AddInEdge(self, src_clust_idx):
		if src_clust_idx not in self.in_edge_list:
			self.in_edge_list.append(src_clust_idx)

	# it adds one in edge information to both the original connectivity (clique) and the final connectivity (tree shape)
	def _AddNoEdge(self, src_clust_idx):
		if src_clust_idx not in self.no_edge_list:
			self.no_edge_list.append(src_clust_idx)
			
	# here the final connectivity is changed (not the original clique based connectivity) -- out edge remove
	def _RemoveOutEdge(self, dest_clust_idx):
		if dest_clust_idx in self.out_edge_list:
			self.out_edge_list.remove(dest_clust_idx)    
		
	# here the final connectivity is changed (not the original clique based connectivity) -- in edge remove
	def _RemoveInEdge(self, dest_clust_idx):
		if dest_clust_idx in self.in_edge_list:
			self.in_edge_list.remove(dest_clust_idx)    
		
	# here the final connectivity is changed (not the original clique based connectivity) -- no edge remove
	def _RemoveNoEdge(self, dest_clust_idx):
		if dest_clust_idx in self.no_edge_list:
			self.no_edge_list.remove(dest_clust_idx)    
		
	def _PrintClusterInfo(self, key, Output_Text_File):
		fp = open(Output_Text_File, 'a')    
		fp.write('\n cluster key: ' + str(key))
		fp.write('\n species list: ' + str(self.Species_List))
		#print 'its indegree: ', self.indegree
		#print 'its outdegree: ', self.outdegree
		fp.write('\n out edge list: ' + str(self.out_edge_list))
		fp.write('\n in edge list: ' + str(self.in_edge_list))
		fp.close()    
    