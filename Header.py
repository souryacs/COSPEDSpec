#!/usr/bin/env python

import dendropy
from dendropy import treesim 
from dendropy import Tree, Taxon, TaxonSet, Node
from dendropy import treecalc
import itertools
import numpy
#from numpy import *
#from numpy import newaxis
import time
import os
import re
from cStringIO import StringIO
from operator import itemgetter
from optparse import OptionParser
from collections import deque
import math
from math import *
#import scipy
#from scipy import optimize 
#from scipy.optimize import minimize

# we define custom edge types
BI_DIRECTED_EDGE = 0	# equality relationship
DIRECTED_OUT_EDGE = 1
DIRECTED_IN_EDGE = 2
NO_EDGE = 3	# no relationship
UNDEFINED_EDGE = 4
ALL_EDGE = 5	# concerns about all edge types

''' this variable is for establish a new connection which is not there in the source trees
that is, two nodes are not in the same source tree, 
thus there is no relationship (even the relationship "NO_EDGE" is not there)
and during tree construction, the connection is coming between them '''
ORIG_DIFF_SRC_TREE_CONNECT_SCORE = 0	#-1

''' this variable signifies that original 'NO_EDGE' connection is not preserved in the derived tree '''
ORIG_NO_EDGE_BECOME_CONNECTED = -3 	

''' this is for the case where original 'NO_EDGE' connection is remained in the derived tree '''
ORIG_NO_EDGE_REMAIN_NO_EDGE = 2


""" declaration of global variables needed to store the structures 
different from the earlier version of supertree program """

""" this is a dictionary storing cluster of nodes 
each cluster is basically a collection of nodes having equality relationship between the nodes """
Cluster_Info_Dict = dict()

""" the dictionary defines one particular taxa 
individual taxa contains the relationship information with other taxa """
Taxa_Info_Dict = dict()

""" this dictionary defines the taxa pair relations
each entry is indexed by two nodes """
TaxaPair_Reln_Dict = dict()

''' this is the list storing the costs of relations for taxa pairs depicting multi relation instance
it stores the taxa pair and the corresponding scores of different relation types
the list is used to extract the maximum score (corresponding to a particular relation) for resolving that taxa pair '''
Cost_List_Taxa_Pair_Multi_Reln = []

''' if SINGLE_EDGE_TYPE_CONN_PRIORITY option is set, then this list stores the cases 
where a pair of taxa is connected only by a single relation type  (with respect to all the candidate source trees)
cost corresponding to that single relation instance is accounted  '''
Cost_List_Taxa_Pair_Single_Reln = [] 

''' this list contains the transitive relationships 
which are established either by original edge connection
or by derived edge connection due to the earlier edge connections
all such edges are maintained in this list
finally, this list is used to update the cost of the remaining non processed edges '''
EDGE_PROCESSED_LIST = []

""" this list contains the complete set of taxa present in the input source trees """
COMPLETE_INPUT_TAXA_LIST = []

""" this list contains the current set of active cluster indices """
CURRENT_CLUST_IDX_LIST = []

# this variable stores the total number of taxa for all the source trees
number_of_taxa = 0

# this variable notes the count of input source trees  
tree_count = 0

# counter of node connection
nodes_connected = 0

# this is the debug level
# set for printing the necessary information
DEBUG_LEVEL = 0

# this text file stores all the printing output
Output_Text_File = 'complete_output_description.txt'

# these two status are for conflict detection routine
# one is when the relation is already established
# other is when the relation is not permitted
RELN_ALREADY_ESTABLISHED = 1
RELN_NOT_PERMITTED = 2

##-----------------------------------------------------
''' this class defines a leaf node of input candidate source trees
that is, it corresponds to one taxa '''
class Single_Taxa(object):
  def __init__(self):
    """ this variable signifies the cluster index that the taxa belongs 
    if the value is -1 (< 0) then it is still not inserted in a cluster 
    otherwise it is part of a valid cluster """
    self.clust_idx_part = -1
    """ these lists store original edge lists (from the input data)
    that is, from the data of input candidate source trees, these lists contain the nodes (taxa) 
    which are related to current taxa 
    lists are classified according to the relationship type
    these lists are updated during the formation of consensus tree, for lowering computation """
    self.orig_out_edge_list = []
    self.orig_in_edge_list = []
    self.orig_eq_edge_list = []
    self.orig_no_edge_list = []
    """ these lists store the final obtained lists (obtained using direct edges) corresponding to the final consensus tree
    these lists are updated according to the consensus tree formation
    lists are classified according to the relationship type """
    self.Final_out_edge_list = []
    self.Final_in_edge_list = []    
    self.Final_eq_edge_list = []
    self.Final_no_edge_list = []

  def _Get_Taxa_Part_Clust_Idx(self):
    return self.clust_idx_part
    
  def _Set_Clust_Idx_taxa_Part(self, inp_clust_idx):
    self.clust_idx_part = inp_clust_idx
  
  # these functions remove taxa from the original edge list of the current taxa
  # according to the original candidate source trees, list of taxa related via certain relationship type were maintained
  # those information is deleted if not required further
  def _SelectivelyRemoveEdgeOrigConn(self, edge_type, other_taxa_label):
    if (edge_type == BI_DIRECTED_EDGE) and (other_taxa_label in self.orig_eq_edge_list):
      self.orig_eq_edge_list.remove(other_taxa_label)
    elif (edge_type == DIRECTED_IN_EDGE) and (other_taxa_label in self.orig_in_edge_list):
      self.orig_in_edge_list.remove(other_taxa_label)
    elif (edge_type == DIRECTED_OUT_EDGE) and (other_taxa_label in self.orig_out_edge_list):
      self.orig_out_edge_list.remove(other_taxa_label)
    elif (edge_type == NO_EDGE) and (other_taxa_label in self.orig_no_edge_list):
      self.orig_no_edge_list.remove(other_taxa_label)
    
  def _RemoveEdgeFromOriginalConnectivity(self, other_taxa_label):
    for edge_type in range(4):
      self._SelectivelyRemoveEdgeOrigConn(edge_type, other_taxa_label)
        
  # these functions add taxa to the original edge list of the current taxa
  # according to the original candidate source trees, list of taxa related via certain relationship type are inserted
  def _AddOrigEdge(self, other_taxa_label, edge_type):
    if (edge_type == BI_DIRECTED_EDGE) and (other_taxa_label not in self.orig_eq_edge_list):
      self.orig_eq_edge_list.append(other_taxa_label)
    if (edge_type == DIRECTED_IN_EDGE) and (other_taxa_label not in self.orig_in_edge_list):
      self.orig_in_edge_list.append(other_taxa_label)
    if (edge_type == DIRECTED_OUT_EDGE) and (other_taxa_label not in self.orig_out_edge_list):
      self.orig_out_edge_list.append(other_taxa_label)
    if (edge_type == NO_EDGE) and (other_taxa_label not in self.orig_no_edge_list):
      self.orig_no_edge_list.append(other_taxa_label)
      
  # this function inserts the final connected taxa with respect to current taxa
  # depending on the connection type, appropriate edge list is considered
  # after inserting the edge type (w.r.t target taxa, those information is removed from the original edge list)
  def _AddFinalEdge(self, edge_type, other_taxa_label):
    flag = 0
    if (edge_type == BI_DIRECTED_EDGE):
      if (other_taxa_label not in self.Final_eq_edge_list):
	self.Final_eq_edge_list.append(other_taxa_label)
	flag = 1
    if (edge_type == DIRECTED_IN_EDGE):
      if (other_taxa_label not in self.Final_in_edge_list):
	self.Final_in_edge_list.append(other_taxa_label)
	flag = 1
    if (edge_type == DIRECTED_OUT_EDGE):
      if (other_taxa_label not in self.Final_out_edge_list):
	self.Final_out_edge_list.append(other_taxa_label)	
	flag = 1
    if (edge_type == NO_EDGE):
      if (other_taxa_label not in self.Final_no_edge_list):
	self.Final_no_edge_list.append(other_taxa_label)	    
	flag = 1
    # now remove this edge information from the original connectivity edges, 
    # as only the remaining non-processed (non-finalized) nodes will be considered
    self._RemoveEdgeFromOriginalConnectivity(other_taxa_label)
    return flag
      
  # these functions return the final connectivity (w.r.t consensus tree) of the current taxa
  # depending on the edge type
  def _GetFinalEqEdgeList(self):
    return self.Final_eq_edge_list
    
  def _GetFinalOutEdgeList(self):
    return self.Final_out_edge_list

  def _GetFinalInEdgeList(self):
    return self.Final_in_edge_list
    
  def _GetFinalNoEdgeList(self):
    return self.Final_no_edge_list    
  
  # these functions return list of taxa w.r.t certain edge types those are yet not considered (or discarded)
  # for final supertree construction
  # depending on the edge type
  def _GetNonProcessedEqEdgeList(self):
    return self.orig_eq_edge_list
    
  def _GetNonProcessedOutEdgeList(self):
    return self.orig_out_edge_list

  def _GetNonProcessedInEdgeList(self):
    return self.orig_in_edge_list

  def _GetNonProcessedNoEdgeList(self):
    return self.orig_no_edge_list
    
  # this function removes a specified taxa from a final edge list
  # selected by the specified edge_type
  def _DelFinalEdge(self, edge_type, other_taxa_label):
    if (edge_type == DIRECTED_OUT_EDGE):
      if other_taxa_label in self.Final_out_edge_list:
	self.Final_out_edge_list.remove(other_taxa_label)
    elif (edge_type == DIRECTED_IN_EDGE):
      if other_taxa_label in self.Final_in_edge_list:
	self.Final_in_edge_list.remove(other_taxa_label)
    elif (edge_type == BI_DIRECTED_EDGE):
      if other_taxa_label in self.Final_eq_edge_list:
	self.Final_eq_edge_list.remove(other_taxa_label)
    else:	#if (edge_type == NO_EDGE):
      if other_taxa_label in self.Final_no_edge_list:
	self.Final_no_edge_list.remove(other_taxa_label)
	
  # this function prints the information of one taxa after processing input candidate source trees
  def _PrintOriginalTaxaInfo(self, key, Output_Text_File):
    fp = open(Output_Text_File, 'a')
    fp.write('\n taxa key: ' + str(key))
    fp.write('\n taxa index in COMPLETE_INPUT_TAXA_LIST: ' + str(COMPLETE_INPUT_TAXA_LIST.index(key)))
    fp.write('\n taxa is part of the cluster ID: ' + str(self.clust_idx_part))
    fp.write('\n original out edge list: ' + str(self.orig_out_edge_list))
    fp.write('\n original in edge list: ' + str(self.orig_in_edge_list))
    fp.write('\n original eq edge list: ' + str(self.orig_eq_edge_list))
    fp.write('\n original no edge list: ' + str(self.orig_no_edge_list))
    fp.close()
    
  # this function is called after formation of consensus tree
  def _PrintFinalTaxaInfo(self, key, Output_Text_File):
    fp = open(Output_Text_File, 'a')    
    fp.write('\n taxa key: ' + str(key))
    fp.write('\n taxa index in COMPLETE_INPUT_TAXA_LIST: ' + str(COMPLETE_INPUT_TAXA_LIST.index(key)))
    fp.write('\n taxa is part of the cluster ID: ' + str(self.clust_idx_part))
    fp.write('\n final out edge list: ' + str(self.Final_out_edge_list))
    fp.write('\n final in edge list: ' + str(self.Final_in_edge_list))
    fp.write('\n final eq edge list: ' + str(self.Final_eq_edge_list))
    fp.write('\n final no edge list: ' + str(self.Final_no_edge_list))
    fp.close()
        
##-----------------------------------------------------
""" this class defines the connectivity relationship between a pair of taxa
initially the information are obtained from the input source trees
later the contents of these class instances are modified according to the generation of the consensus tree
key of this class --- taxa1, taxa2  
in the class, the edge type signifies the relation between a pair of taxa
"""
class Reln_TaxaPair(object):
  def __init__(self, fin_edge_tp):    
    # this is the final selected edge type that is established between two taxa 
    self.final_selected_edge_type = fin_edge_tp
    ''' this variable denotes the no of occurrences of a particular edge type 
    there are 4 types of edges (relationship) between a pair of taxa '''
    self.edge_weight = [0] * 4    
    ''' a connection priority value is defined as the 
    no of occurrences of this particular edge type between these pair of nodes 
    minus the sum of no of occurrences of other edge types between these pair of nodes '''
    self.conn_pr_val = [0] * 4    
    ''' this cost variable denotes the cost associated with different types of edge connection between 
    these pair of nodes considered 
    this value is updated during generation of the consensus tree '''
    self.Connect_Edge_Cost = [0] * 4
    ''' this variable stores the number of levels that belong within this taxa pair
    it is the sum of branches existing between individual nodes to their MRCA 
    for all the trees, the sum is accumulated '''
    self.level_sum_all_trees = 0
    ''' this variable stores the no of trees supporting the taxa pair '''
    self.supporting_trees = 0
    
  def _AddLevels(self, incr_level):
    self.level_sum_all_trees = self.level_sum_all_trees + incr_level
    
  def _AddSupportingTree(self):
    self.supporting_trees = self.supporting_trees + 1
    
  # modified - sourya
  def _GetAvgLevelPerTree(self):
    return (self.level_sum_all_trees / self.supporting_trees)
    
  def _GetEdgeWeight(self, edge_type):
    return self.edge_weight[edge_type]      
    
  def _GetEdgeCost_ConnReln(self, edge_type):
    return self.Connect_Edge_Cost[edge_type]
    
  def _IncrEdgeCost_ConnReln(self, edge_type, incr_cost):
    self.Connect_Edge_Cost[edge_type] = self.Connect_Edge_Cost[edge_type] + incr_cost

  # this function adds one edge count (with a given input edge type)
  def _AddEdgeCount(self, edge_type, no_of_levels = 0):
    # increment the frequency of particular relation type between this taxa pair
    if (no_of_levels == 0):
      self.edge_weight[edge_type] = self.edge_weight[edge_type] + 1
    else:
      self.edge_weight[edge_type] = self.edge_weight[edge_type] + no_of_levels
    
  # this function appends one final connected edge (for use in final supertree) 
  def _UpdateFinalEdgeInfo(self, edge_type):
    self.final_selected_edge_type = edge_type
    
  # this function prints the relationship information
  def _PrintRelnInfo(self, key, Output_Text_File):
    fp = open(Output_Text_File, 'a')    
    fp.write('\n taxa pair key: ' + str(key))
    fp.write('\n edges [type/count/conn_pr_val/score]: ')
    for i in range(4):
      fp.write('\n [' + str(i) + '/' + str(self.edge_weight[i]) + '/' + str(self.conn_pr_val[i]) + '/' + str(self.Connect_Edge_Cost[i]) + ']')
    fp.write('\n final selected edge type : ' + str(self.final_selected_edge_type))
    fp.close()
         
  # this function computes the score metric value associated with individual pair of taxa 
  def _SetCostMetric(self):      
    for edge_type in range(4):
      # assign the score metric for this edge type
      self.Connect_Edge_Cost[edge_type] = self.edge_weight[edge_type] * self.conn_pr_val[edge_type]
      
  # this function returns the connection priority value for input edge type
  def _GetConnPrVal(self, edge_type):
    return self.conn_pr_val[edge_type]
      
  ''' this function calculates connection priority value for each of the edge types, 
  for this particular connection between a pair of nodes in the final tree '''
  def _SetConnPrVal(self, single_edge_prior):
    # this is the sum of all the edge type instances (no of occurrences)
    listsum = sum(self.edge_weight)
    # now determine the connection priority of a particular edge type with respect to other edges     
    for edge_type in range(4):
      # here we use the difference of current edge type frequency with the frequencies of all other edge types 
      self.conn_pr_val[edge_type] = 2 * self.edge_weight[edge_type] - listsum
    ''' this code section is used when there exists NO EDGE relationship between a pair of taxa
    and we want to detect it '''
    if (not single_edge_prior):
      """ if there is no vote for any particular edge type other than NO_EDGE,
      (that is, corresponding settings did not occur in any of the source tree)
      then we make only the NO_EDGE settings as valid - 
      they will only be considered for joining this pair in the final tree """
      if (self.edge_weight[NO_EDGE] != 0)\
	and (self.edge_weight[DIRECTED_IN_EDGE] == 0)\
	and (self.edge_weight[DIRECTED_OUT_EDGE] == 0)\
	and (self.edge_weight[BI_DIRECTED_EDGE] == 0):
	return 1
      else:
	return 0
    else:
      outlist = [0, NO_EDGE]
      for edge_type in range(4):
	if (self.edge_weight[edge_type] == listsum) and (listsum > 0):
	  outlist = [1, edge_type]
	  break
	elif (self.edge_weight[edge_type] > 0) and (self.edge_weight[edge_type] < listsum):
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
    