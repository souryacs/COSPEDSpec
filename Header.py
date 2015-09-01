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
BI_DIRECTED_EDGE = 0	# relation r3
DIRECTED_OUT_EDGE = 1	# relation r1
DIRECTED_IN_EDGE = 2	# relation r2
NO_EDGE = 3	# relation r4
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
    self.edge_weight = [0] * 4    
    
    """ 
    a connection priority value is defined as the 
    no of occurrences of this particular relation between this pair of taxa 
    minus the sum of no of occurrences of other relation types between this couplet
    """
    self.conn_pr_val = [0] * 4    
    
    """ 
    this is the support score for different types of relations between a couplet
    """
    self.Connect_Edge_Cost = [0] * 4
    
    """ 
    this variable stores the no of trees supporting the taxa pair 
    """
    self.supporting_trees = 0
    
    """ 
    For this couplet, it stores the extra gene count with respect to all the gene trees
    """
    self.XL_sum_gene_trees = 0
    
  def _AddXLVal(self, XL_val):
    self.XL_sum_gene_trees = self.XL_sum_gene_trees + XL_val
    
  def _GetXLSumGeneTrees(self):
    return self.XL_sum_gene_trees
  
  def _GetNormalizedXLSumGeneTrees(self):
    return (self.XL_sum_gene_trees * 1.0) / self.supporting_trees
    
  def _AddSupportingTree(self):
    self.supporting_trees = self.supporting_trees + 1
        
  def _GetEdgeWeight(self, edge_type):
    return self.edge_weight[edge_type]      
    
  def _GetEdgeCost_ConnReln(self, edge_type):
    return self.Connect_Edge_Cost[edge_type]
    
  def _IncrEdgeCost_ConnReln(self, edge_type, incr_cost):
    self.Connect_Edge_Cost[edge_type] = self.Connect_Edge_Cost[edge_type] + incr_cost

  # this function adds one frequency count (with a given input relation type)
  def _AddEdgeCount(self, edge_type):
    self.edge_weight[edge_type] = self.edge_weight[edge_type] + 1
    
  # this function prints the relationship information
  def _PrintRelnInfo(self, key, Output_Text_File):
    fp = open(Output_Text_File, 'a')    
    fp.write('\n taxa pair key: ' + str(key))
    fp.write('\n edges [type/count/conn_pr_val/score]: ')
    for i in range(4):
      fp.write('\n [' + str(i) + '/' + str(self.edge_weight[i]) + '/' + str(self.conn_pr_val[i]) + '/' + str(self.Connect_Edge_Cost[i]) + ']')
    fp.write('\n Sum of extra lineage : ' + str(self.XL_sum_gene_trees))
    fp.write('\n No of supporting trees : ' + str(self.supporting_trees))
    fp.write('\n Normalized XL sum : ' + str(self._GetNormalizedXLSumGeneTrees()))
    fp.close()
         
  # this function computes the support score metric value associated with individual pair of taxa 
  def _SetCostMetric(self):      
    for edge_type in range(4):
      # assign the score metric for this edge type
      self.Connect_Edge_Cost[edge_type] = self.edge_weight[edge_type] * self.conn_pr_val[edge_type]
      
  # this function returns the connection priority value for input relation 
  def _GetConnPrVal(self, edge_type):
    return self.conn_pr_val[edge_type]
      
  """ 
  this function calculates connection priority value for each of the relation types, 
  """
  def _SetConnPrVal(self, single_edge_prior):
    # this is the sum of frequencies for all the relation types
    listsum = sum(self.edge_weight)
    # now determine the connection priority of a particular relation type with respect to other relations     
    for edge_type in range(4):
      # here we use the difference of current relation type frequency with the frequencies of all other relations
      self.conn_pr_val[edge_type] = 2 * self.edge_weight[edge_type] - listsum
    """ 
    this code section is used when there exists an unique relation (non-conflicting couplet)
    and we try to detect it
    """
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
    