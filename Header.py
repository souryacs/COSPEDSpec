#!/usr/bin/env python

import dendropy
from dendropy import Tree, Taxon, TaxonSet, Node
import numpy
import time
import os
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
this list contains the complete set of taxa present in the input trees 
"""
COMPLETE_INPUT_TAXA_LIST = []

""" 
this list contains the current set of active taxa cluster (indices)
"""
CURRENT_CLUST_IDX_LIST = []

# this is the debug level
# set for printing the necessary information
DEBUG_LEVEL = 2

MAJORITY_CONSENSUS_RATIO = 0.6
LEVEL_COUNT_VAL_CONSENSUS_RATIO = 0.7

# variables used to denote whether we use traditional NJ method
# or use a variant of it, namely the agglomerative clustering
TRADITIONAL_NJ = 1
AGGLO_CLUST = 2

# add - sourya
MODE_PERCENT = 0.4	#0.5	#0.35
MODE_BIN_COUNT = 40

"""
this is a threshold corresponding to the selection of R1 or R2 relation
for a non conflicting couplet with negative support score of the corresponding relation
"""
#R1R2Reln_MAJ_THRS = 0.7
R1R2Reln_MAJ_THRS_high = 0.75	#0.8
R1R2Reln_MAJ_THRS_low = 0.65	#0.7	#0.8
R1R2Reln_MAJ_THRS_very_low = 0.6	#0.65

"""
this is a threshold corresponding to the selection of R3 relation
for a non conflicting couplet with negative support score of the corresponding relation
"""
R3Reln_MAJ_THRS = 0.2
#R3Reln_MAJ_THRS_LOW = 0.15

"""
for a couplet, relations with frequency of this percentage of the max (consensus) frequency
will be included in the score queue
"""
# it should be fixed at 0.35 - sourya
PERCENT_MAX_FREQ = 0.35	#0.5	#sourya

#CONSENSUS_FREQ_RATIO_THR = 0.7	#0.8

#PRIORITY_PERCENT = 0.1

"""
this list contains the set of clusters 
which need to be processed for setting at least one directed out edge to another cluster
"""
Candidate_Out_Edge_Cluster_List = []

"""
this threshold is applied at the last stage
when a sibling is broken in R1 / R2 relations 
based on the XL value distribution
"""
#VARIATION_XL_THR = 0.08	#0.25

#-----------------------------------------------------
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
		# modified - sourya
		self.XL_sum_gene_trees = []	#0
		
		"""
		for a couplet xy, and for a relation r3, following array has 3 elements:
		1) count when level of x < level of y (relation r1 similar)
		2) count when level of x > level of y (relation r2 similar)
		3) count when level of x = level of y (relation r3 similar)
		"""
		self.ALL_Reln_Level_Diff_Info_Count = [0] * 3
		self.ALL_Reln_Level_Diff_Val_Count = [0] * 3
		
		"""
		this is a variable containing the binned average of the XL values
		of very high frequency
		initially the value is set as -1, to signify that the computation is not done
		once the computation (for a couplet) is done, the value is subsequently used and returned
		"""
		self.binned_avg_XL = -1
		
		"""
		this is a list containing the number of instances
		when R4 relation is actually a pseudo R1 relation
		"""
		self.freq_R4_pseudo_R1R2 = [0] * 2
		"""
		this is a list of allowed relations between this couplet
		we allow only those relations which have a significant frequency compared to the number of supporting trees
		"""
		self.allowed_reln_list = []
	
	#----------------------------------
	"""
	this function returns the ratio of level value
	@param: idx: 	if 0, returns ratio w.r.t r1
							  if 1, returns ratio w.r.t r2
	"""
	def _GetLevelValRatio(self, idx):
		level_val_r1 = self.ALL_Reln_Level_Diff_Val_Count[0]
		level_val_r2 = self.ALL_Reln_Level_Diff_Val_Count[1]
		if (idx == 0):
			if ((level_val_r1 + level_val_r2) > 0):
				return (level_val_r1 * 1.0) / (level_val_r1 + level_val_r2)
			else:
				return 0
		else:	#if (idx == 1):
			if ((level_val_r1 + level_val_r2) > 0):
				return (level_val_r2 * 1.0) / (level_val_r1 + level_val_r2)
			else:
				return 0
			
		return 0
	
	"""
	this function checks whether for a conflicting taxa pair with negative support score
	R3 relation can be applied between this couplet
	"""
	def _Check_Reln_R3_Majority(self, outfile):
		level_val_r1 = self.ALL_Reln_Level_Diff_Val_Count[0]
		level_val_r2 = self.ALL_Reln_Level_Diff_Val_Count[1]
		lev_diff = math.fabs((level_val_r2 - level_val_r1) * 1.0)
		fr3 = self.freq_count[RELATION_R3]
		
		if (DEBUG_LEVEL >= 2):
			fp = open(outfile, 'a')
			fp.write('\n fr3: ' + str(fr3) + ' supporting trees : ' + str(self.supporting_trees))
			if ((level_val_r2 + level_val_r1) > 0):
				fp.write(' Ratio: ' + str(round((lev_diff / (level_val_r2 + level_val_r1)), 2)))  
			fp.close()
		
		"""
		if R3 relation is consensus then we impose R3 relation between this couplet
		"""
		if (self._CheckTargetRelnConsensus(RELATION_R3) == True):
			if (DEBUG_LEVEL >= 2):
				fp = open(outfile, 'a')
				fp.write('\n *** R3 RELATI0N IS THE MAJORITY *** ')  
				fp.close()
			return True
		
		# add - sourya
		"""
		if the excess gene count statistic (mode) is 0 for this couplet
		then surely R3 relation is predominant
		"""
		if (self._GetMultiModeXLVal() == 0):
			if (DEBUG_LEVEL >= 2):
				fp = open(outfile, 'a')
				fp.write('\n *** R3 RELATI0N IS THE MAJORITY *** ')  
				fp.close()
			return True
		# end add - sourya
		
		## add - sourya
		#if ((level_val_r2 + level_val_r1) > 0):
			#if ((round((lev_diff / (level_val_r2 + level_val_r1)), 2)) <= R3Reln_MAJ_THRS_LOW):
				## the level difference should be very small
				#if (DEBUG_LEVEL >= 2):
					#fp = open(outfile, 'a')
					#fp.write('\n *** R3 RELATI0N IS THE MAJORITY *** ')  
					#fp.close()
				#return True
			#if ((round((lev_diff / (level_val_r2 + level_val_r1)), 2)) <= R3Reln_MAJ_THRS) and (fr3 > min(self.freq_count)):
				## the level difference should be very small
				#if (DEBUG_LEVEL >= 2):
					#fp = open(outfile, 'a')
					#fp.write('\n *** R3 RELATI0N IS THE MAJORITY *** ')  
					#fp.close()
				#return True
		#else:	#if ((level_val_r2 + level_val_r1) == 0):
			#if (DEBUG_LEVEL >= 2):
				#fp = open(outfile, 'a')
				#fp.write('\n *** R3 RELATI0N IS THE MAJORITY *** ')  
				#fp.close()
			#return True
		## end add - sourya
		
		## comment - sourya
		"""
		R3 relation should not be the relation with minimum frequency among all constituent relations
		unless the frequency of R3 relation is at least 20% of the input supported trees
		this 20\% threshold is empirical
		"""
		# modified - sourya
		#if (fr3 > min(self.freq_count)) or (fr3 >= (0.2 * self.supporting_trees)):
		#if (fr3 > min(self.freq_count)) and (fr3 >= (0.6 * max(self.freq_count))):
		if (fr3 > min(self.freq_count)):
			if ((level_val_r2 + level_val_r1) > 0):
				if ((round((lev_diff / (level_val_r2 + level_val_r1)), 2)) <= R3Reln_MAJ_THRS):
					# the level difference should be very small
					if (DEBUG_LEVEL >= 2):
						fp = open(outfile, 'a')
						fp.write('\n *** R3 RELATI0N IS THE MAJORITY *** ')  
						fp.close()
					return True
			else:	#if ((level_val_r2 + level_val_r1) == 0):
				if (DEBUG_LEVEL >= 2):
					fp = open(outfile, 'a')
					fp.write('\n *** R3 RELATI0N IS THE MAJORITY *** ')  
					fp.close()
				return True
		## end comment - sourya

		if (DEBUG_LEVEL >= 2):
			fp = open(outfile, 'a')
			fp.write('\n R3 RELATI0N IS NOT THE MAJORITY')  
			fp.close()
		
		return False
		
	#----------------------------------
	"""
	these functions adjust the set of allowed relations among a couplet
	"""
	def _AddAllowedReln(self, inp_reln):
		self.allowed_reln_list.append(inp_reln)
		
	def _GetAllowedRelnList(self):
		return self.allowed_reln_list
	
	def _RemoveAllowedReln(self, inp_reln):
		if inp_reln in self.allowed_reln_list:
			self.allowed_reln_list.remove(inp_reln)
	#----------------------------------	
	
	def _AddFreqPseudoR1(self, idx, r=1):
		self.freq_R4_pseudo_R1R2[idx] = self.freq_R4_pseudo_R1R2[idx] + r
		
	def _GetFreqPseudoR1(self, idx):
		return self.freq_R4_pseudo_R1R2[idx]
		
	def _NormalizeR1R2LevelDiff(self):
		for i in range(3):
			self.ALL_Reln_Level_Diff_Val_Count[i] = (self.ALL_Reln_Level_Diff_Val_Count[i] * 1.0) / self.supporting_trees
		
	"""
	this function computes the level difference (relation based) count 
	# parameters:
	@src_reln: corresponds to the level of relation from which the subtraction will take place
	@dest_reln: if between 0 and 2, corresponds to the subtracted relation
	if it is -1, all relations except the src_reln will be subtracted
	@abs_comp: if true, will compute the absolute value of the subtracted quantity
	@norm: if true, will normalize the subtracted value with the number of supporting trees
	"""
	def _GetRelnLevelDiff(self, src_reln, dest_reln, abs_comp, norm):
		reln_array = [RELATION_R1, RELATION_R2, RELATION_R3]
		src_reln_idx = reln_array.index(src_reln)

		target_val = self.ALL_Reln_Level_Diff_Info_Count[src_reln_idx]
		
		if (dest_reln == -1):
			# all relations will be subtracted
			for dest_reln_idx in range(3):
				if (dest_reln_idx == src_reln_idx):
					continue
				target_val = target_val - self.ALL_Reln_Level_Diff_Info_Count[dest_reln_idx]
		else:
			dest_reln_idx = reln_array.index(dest_reln)
			target_val = target_val - self.ALL_Reln_Level_Diff_Info_Count[dest_reln_idx]
			
		if (abs_comp == True):
			target_val = math.fabs(target_val)
			
		if (norm == True):
			target_val = (target_val * 1.0) / self.supporting_trees
			
		return target_val

	"""
	this function checks whether the input relation type is a consensus relation
	among this couplet
	"""
	def _CheckTargetRelnConsensus(self, inp_reln_type):
		if (self.freq_count[inp_reln_type] == max(self.freq_count)):
			return True
		return False
	
	def _CheckTargetRelnLevelConsensus(self, src_reln, only_cons=0):
		reln_array = [RELATION_R1, RELATION_R2, RELATION_R3]
		src_reln_idx = reln_array.index(src_reln)
		sum_level_count = sum(self.ALL_Reln_Level_Diff_Info_Count)
		sum_level_val = self.ALL_Reln_Level_Diff_Val_Count[0] + self.ALL_Reln_Level_Diff_Val_Count[1]
		if (src_reln_idx == 0) or (src_reln_idx == 1):
			if (only_cons == 0):
				if (self.ALL_Reln_Level_Diff_Info_Count[src_reln_idx] >= (MAJORITY_CONSENSUS_RATIO * sum_level_count)) and \
					(self.ALL_Reln_Level_Diff_Val_Count[src_reln_idx] >= (LEVEL_COUNT_VAL_CONSENSUS_RATIO * sum_level_val)):
					return 1
			else:
				if (self.ALL_Reln_Level_Diff_Info_Count[src_reln_idx] > (0.5 * sum_level_count)) and \
					(self.ALL_Reln_Level_Diff_Val_Count[src_reln_idx] > (0.5 * sum_level_val)):
					return 1
		else:
			if (self.ALL_Reln_Level_Diff_Info_Count[2] > (self.ALL_Reln_Level_Diff_Info_Count[0] + self.ALL_Reln_Level_Diff_Info_Count[1])):
				return 1
			
		return 0
		
	def _IncrAllRelnLevelDiffInfoCount(self, idx, val):
		self.ALL_Reln_Level_Diff_Info_Count[idx] = self.ALL_Reln_Level_Diff_Info_Count[idx] + 1
		self.ALL_Reln_Level_Diff_Val_Count[idx] = self.ALL_Reln_Level_Diff_Val_Count[idx] + val
		
	#------------------------------------------
	def _GetXLList(self):
		return self.XL_sum_gene_trees
	
	"""
	this function adds one XL value, computed for a particular input tree
	to the list of XL values for this couplet
	"""
	def _AddXLVal(self, XL_val):
		#self.XL_sum_gene_trees = self.XL_sum_gene_trees + XL_val
		self.XL_sum_gene_trees.append(XL_val)	# modified - sourya
		
	"""
	this function returns the list of XL values for this couplet
	with respect to individual input trees
	"""
	def _GetXLSumGeneTrees(self):
		#return self.XL_sum_gene_trees
		return sum(self.XL_sum_gene_trees)	# modified - sourya
	
	"""
	this function computes the average of XL measures
	"""
	def _GetAvgXLGeneTrees(self):
		#return (self.XL_sum_gene_trees * 1.0) / self.supporting_trees
		return (sum(self.XL_sum_gene_trees) * 1.0) / self.supporting_trees	# modified - sourya
	
	"""
	this function computes the median of XL measures
	"""
	def _GetMedianXLGeneTrees(self):
		return numpy.median(numpy.array(self.XL_sum_gene_trees)) # modified - sourya

	def _GetStdevXL(self):
		return numpy.std(numpy.array(self.XL_sum_gene_trees))
	
	def _GetVarianceXL(self):
		return numpy.var(numpy.array(self.XL_sum_gene_trees))
	
	"""
	function to return the average of XL values for this couplet
	depending on the user parameters, average, median, or binned average XL is returned
	"""
	def _GetNormalizedXLSumGeneTrees(self, dist_type):
		if (dist_type == 1):
			return self._GetAvgXLGeneTrees()
		elif (dist_type == 2):
			return self._GetMedianXLGeneTrees()
		elif (dist_type == 3):
			return self._GetMultiModeXLVal()
		elif (dist_type == 4):
			return min(self._GetAvgXLGeneTrees(), self._GetMedianXLGeneTrees())
		elif (dist_type == 5):
			# modified - sourya
			#return min(self._GetAvgXLGeneTrees(), self._GetMedianXLGeneTrees(), self._GetMultiModeXLVal())
			# we return average of these three quantities
			return (self._GetAvgXLGeneTrees() + self._GetMedianXLGeneTrees() + self._GetMultiModeXLVal()) / 3.0
		elif (dist_type == 6):
			return min(self._GetMedianXLGeneTrees(), self._GetMultiModeXLVal())
		
	"""
	this function computes the binned average of XL values associated for this couplet
	"""
	def _GetMultiModeXLVal(self, Output_Text_File=None):
		if (self.binned_avg_XL == -1):
			
			Bin_Width = (1.0 / MODE_BIN_COUNT)
			len_list = [0] * MODE_BIN_COUNT
			
			if Output_Text_File is not None:
				fp = open(Output_Text_File, 'a') 
			
			# sort the XL list
			self.XL_sum_gene_trees.sort()
			
			for j in range(len(self.XL_sum_gene_trees)):
				curr_xl_val = self.XL_sum_gene_trees[j]
				bin_idx = int(curr_xl_val / Bin_Width)
				if (bin_idx == MODE_BIN_COUNT):
					bin_idx = bin_idx - 1
				len_list[bin_idx] = len_list[bin_idx] + 1
			
			if Output_Text_File is not None:
				for i in range(MODE_BIN_COUNT):
					fp.write('\n bin idx: ' + str(i) + ' len:  ' + str(len_list[i]))
			
			# this is the maximum length of a particular bin
			# corresponding to max frequency
			max_freq = max(len_list)
			
			if Output_Text_File is not None:
				fp.write('\n Max freq: ' + str(max_freq))
			
			num = 0
			denom = 0
			for i in range(MODE_BIN_COUNT):
				if (len_list[i] >= (MODE_PERCENT * max_freq)):
					list_start_idx = sum(len_list[:i])
					list_end_idx = list_start_idx + len_list[i] - 1
					value_sum = sum(self.XL_sum_gene_trees[list_start_idx:(list_end_idx+1)])
					num = num + value_sum
					denom = denom + len_list[i]
					if Output_Text_File is not None:
						fp.write('\n Included bin idx: ' + str(i) + ' starting point: ' + str(list_start_idx) \
							+ 'ending point: ' + str(list_end_idx) + ' sum: ' + str(value_sum))
			
			self.binned_avg_XL = (num / denom)
			
			if Output_Text_File is not None:
				fp.write('\n Final binned average XL: ' + str(self.binned_avg_XL))
				fp.close()
			
		return self.binned_avg_XL
	
	#------------------------------------------
	"""
	this function adds one supporting tree for this couplet
	"""
	def _AddSupportingTree(self):
		self.supporting_trees = self.supporting_trees + 1
	
	"""
	this function returns the number of input trees supporting this couplet
	"""
	def _GetNoSupportTrees(self):
		return self.supporting_trees
	#------------------------------------------
	"""
	this function returns the frequency of the consensus relation
	"""
	def _GetConsensusFreq(self):
		return max(self.freq_count)
	
	"""
	this function returns the frequency of the input relation
	specified by the variable 'reln_type'
	"""
	def _GetEdgeWeight(self, reln_type):
		return self.freq_count[reln_type]      

	"""
	this function adds a specified frequency count (default 1)
	corresponding to the relation specified by the variable 'reln_type'
	"""
	def _AddEdgeCount(self, reln_type, r=1):
		# modified - sourya
		if 0:	#(reln_type == RELATION_R3):
			self.freq_count[reln_type] = self.freq_count[reln_type] + (2 * r)
		else:
			self.freq_count[reln_type] = self.freq_count[reln_type] + r
	#------------------------------------------
	"""
	this function returns the support score of the input relation
	specified by the variable 'reln_type'
	"""
	def _GetEdgeCost_ConnReln(self, reln_type):
		return self.support_score[reln_type]

	""" 
	this function computes the support score value associated with individual couplet
	for all different relations
	"""
	def _SetCostMetric(self):      
		for reln_type in range(4):
			# assign the score metric for this edge type
			self.support_score[reln_type] = self.freq_count[reln_type] * self.priority_reln[reln_type]

	"""
	this function updates (increments) the support score of the input relation
	specified by the variable 'reln_type'
	and by the amount 'incr_cost'
	"""
	def _IncrEdgeCost_ConnReln(self, reln_type, incr_cost):
		self.support_score[reln_type] = self.support_score[reln_type] + incr_cost
	#------------------------------------------
	"""
	this function prints all the information associated with a couplet
	"""
	def _PrintRelnInfo(self, key, Output_Text_File=None):
		
		if Output_Text_File is not None:
		
			fp = open(Output_Text_File, 'a')    
			fp.write('\n\n\n taxa pair key: ' + str(key))
			fp.write('\n relations [type/count/priority_reln/score]: ')
			for i in range(4):
				fp.write('\n [' + str(i) + '/' + str(self.freq_count[i]) + '/' + str(self.priority_reln[i]) + '/' + str(self.support_score[i]) + ']')
			#fp.write('\n Sum of extra lineage **** : ' + str(self.XL_sum_gene_trees))
			# add - sourya
			#fp.write('\n XL list **** : ' + str(self.XL_sum_gene_trees))
			fp.write('\n AVERAGE Sum of extra lineage **** : ' + str(self._GetAvgXLGeneTrees()))
			fp.write('\n MEDIAN Sum of extra lineage **** : ' + str(self._GetMedianXLGeneTrees()))
			fp.write('\n Mode Sum of extra lineage **** : ' + str(self._GetMultiModeXLVal()))
			fp.write('\n Mean(Avg + Mode) of extra lineage **** : ' + str((self._GetMultiModeXLVal() + self._GetAvgXLGeneTrees())/2.0))
			fp.write('\n Variance of excess lineage: ' + str(self._GetVarianceXL()))
			fp.write('\n Standard deviation of excess lineage: ' + str(self._GetStdevXL()))
			# end add - sourya
			fp.write('\n No of supporting trees : ' + str(self.supporting_trees))
			#fp.write('\n Normalized XL sum : ' + str(self._GetNormalizedXLSumGeneTrees()))
			fp.write('\n ALL relation based Level diff info count (r1/r2/r3): ' + str(self.ALL_Reln_Level_Diff_Info_Count))
			fp.write('\n ALL relation based Level diff Val count (r1/r2/r3): ' + str(self.ALL_Reln_Level_Diff_Val_Count))
			fp.write('\n R4 relation pseudo (R1/R2) count: ' + str(self.freq_R4_pseudo_R1R2))
			fp.write('\n Allowed relation list: ' + str(self.allowed_reln_list))
			if ((self.ALL_Reln_Level_Diff_Val_Count[0] + self.ALL_Reln_Level_Diff_Val_Count[1]) > 0):
				fp.write('\n Level Val r1 reln ratio : ' + str((self.ALL_Reln_Level_Diff_Val_Count[0] * 1.0) / (self.ALL_Reln_Level_Diff_Val_Count[0] + self.ALL_Reln_Level_Diff_Val_Count[1])))
				fp.write('\n Level Val r2 reln ratio : ' + str((self.ALL_Reln_Level_Diff_Val_Count[1] * 1.0) / (self.ALL_Reln_Level_Diff_Val_Count[0] + self.ALL_Reln_Level_Diff_Val_Count[1])))
				fp.write('\n Level Val r3/r4 reln ratio : ' + str((math.fabs(self.ALL_Reln_Level_Diff_Val_Count[0] - self.ALL_Reln_Level_Diff_Val_Count[1])) / (self.ALL_Reln_Level_Diff_Val_Count[0] + self.ALL_Reln_Level_Diff_Val_Count[1])))
			#if (sum(self.ALL_Reln_Level_Diff_Info_Count) > 0):
				#fp.write('\n Level Count r1 reln ratio : ' + str((self.ALL_Reln_Level_Diff_Info_Count[0] * 1.0) / (sum(self.ALL_Reln_Level_Diff_Info_Count))))
				#fp.write('\n Level Count r2 reln ratio : ' + str((self.ALL_Reln_Level_Diff_Info_Count[1] * 1.0) / (sum(self.ALL_Reln_Level_Diff_Info_Count))))
				#fp.write('\n Level Count r3 reln ratio : ' + str((self.ALL_Reln_Level_Diff_Info_Count[2] * 1.0) / (sum(self.ALL_Reln_Level_Diff_Info_Count))))

			fp.close()
			
		else:
			
			sys.stdout.write('\n\n\n taxa pair key: ' + str(key))
			sys.stdout.write('\n relations [type/count/priority_reln/score]: ')
			for i in range(4):
				sys.stdout.write('\n [' + str(i) + '/' + str(self.freq_count[i]) + '/' + str(self.priority_reln[i]) + '/' + str(self.support_score[i]) + ']')
			sys.stdout.write('\n AVERAGE Sum of extra lineage **** : ' + str(self._GetAvgXLGeneTrees()))
			sys.stdout.write('\n MEDIAN Sum of extra lineage **** : ' + str(self._GetMedianXLGeneTrees()))
			sys.stdout.write('\n Mode Sum of extra lineage **** : ' + str(self._GetMultiModeXLVal()))
			sys.stdout.write('\n No of supporting trees : ' + str(self.supporting_trees))
			sys.stdout.write('\n ALL relation based Level diff info count (r1/r2/r3): ' + str(self.ALL_Reln_Level_Diff_Info_Count))
			sys.stdout.write('\n ALL relation based Level diff Val count (r1/r2/r3): ' + str(self.ALL_Reln_Level_Diff_Val_Count))
			sys.stdout.write('\n R4 relation pseudo (R1/R2) count: ' + str(self.freq_R4_pseudo_R1R2))
			sys.stdout.write('\n Allowed relation list: ' + str(self.allowed_reln_list))
			if ((self.ALL_Reln_Level_Diff_Val_Count[0] + self.ALL_Reln_Level_Diff_Val_Count[1]) > 0):
				sys.stdout.write('\n Level Val r1 reln ratio : ' + str((self.ALL_Reln_Level_Diff_Val_Count[0] * 1.0) / (self.ALL_Reln_Level_Diff_Val_Count[0] + self.ALL_Reln_Level_Diff_Val_Count[1])))
				sys.stdout.write('\n Level Val r2 reln ratio : ' + str((self.ALL_Reln_Level_Diff_Val_Count[1] * 1.0) / (self.ALL_Reln_Level_Diff_Val_Count[0] + self.ALL_Reln_Level_Diff_Val_Count[1])))
				sys.stdout.write('\n Level Val r3/r4 reln ratio : ' + str((math.fabs(self.ALL_Reln_Level_Diff_Val_Count[0] - self.ALL_Reln_Level_Diff_Val_Count[1])) / (self.ALL_Reln_Level_Diff_Val_Count[0] + self.ALL_Reln_Level_Diff_Val_Count[1])))
			sys.stdout.write('\n\n\n')

	#------------------------------------------
	"""
	this function returns the priority value for a given input relation 
	specified by the variable 'reln_type'
	"""
	def _GetConnPrVal(self, reln_type):
		return self.priority_reln[reln_type]
	
	""" 
	this function calculates connection priority value for each of the relation types
	"""
	def _SetConnPrVal(self):
		# this is the sum of frequencies for all the relation types
		listsum = sum(self.freq_count)
		# now determine the connection priority of a particular relation type with respect to other relations     
		for reln_type in range(4):
			# here we use the difference of current relation type frequency with the frequencies of all other relations
			self.priority_reln[reln_type] = 2 * self.freq_count[reln_type] - listsum
		
	#------------------------------------------
	"""
	this function checks whether a couplet is non-conflicting
	that is, only one relation between them exists throughout all the gene trees
	in such a case, a binary variable 1 and the corresponding relation type is returned in the form of a list
	otherwise, a value 0 and a defaukt relation R4 is returned
	"""
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
		# this is a possible R1 reln list
		self.possible_R1_list = []
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
	
	def _GetCardinality(self):
		return len(self.Species_List)
				
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
			
	#--------------------------------------------------------
	# add - sourya
	def _AddPossibleR1(self, dest_clust_idx):
		if dest_clust_idx not in self.possible_R1_list:
			self.possible_R1_list.append(dest_clust_idx)
	
	def _RemovePossibleR1(self, dest_clust_idx):
		if dest_clust_idx in self.possible_R1_list:
			self.possible_R1_list.remove(dest_clust_idx)
			
	def _GetPossibleR1List(self):
		return self.possible_R1_list
	#--------------------------------------------------------
		
	def _PrintClusterInfo(self, key, Output_Text_File):
		fp = open(Output_Text_File, 'a')    
		fp.write('\n cluster key: ' + str(key))
		fp.write('\n species list: ' + str(self.Species_List))
		#print 'its indegree: ', self.indegree
		#print 'its outdegree: ', self.outdegree
		fp.write('\n out edge list: ' + str(self.out_edge_list))
		fp.write('\n in edge list: ' + str(self.in_edge_list))
		fp.write('\n Possible R1 list: ' + str(self.possible_R1_list))	# add - sourya
		fp.close()    
    