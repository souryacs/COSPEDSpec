#!/usr/bin/env python

import Header
from Header import *
import UtilFunc
from UtilFunc import *


##-----------------------------------------------------
# this function updates the list containing the score of node pairs, by an amount delta_score
def UpdateEdgeCostInMemList(Cost_List_Node_Pair, nodeA, nodeB, delta_score, connect_edge_type):
  if (delta_score != 0):
    # form the list key to search the content of the Cost_List_Node_Pair
    key = (nodeA, nodeB)
    edge_cost_orig = TaxaPair_Reln_Dict[key]._GetEdgeCost_ConnReln(connect_edge_type)
    orig_list_key = [nodeA, nodeB, connect_edge_type, edge_cost_orig]
    
    # this condition is to prevent key error 
    # some keys may not be found due to late addition of unnecessary edges (connections)
    if orig_list_key in Cost_List_Node_Pair:
      if (delta_score > 0):
	Queue_Increase_Score(Cost_List_Node_Pair, Cost_List_Node_Pair.index(orig_list_key), (edge_cost_orig + delta_score))
      else:
	Queue_Decrease_Score(Cost_List_Node_Pair, Cost_List_Node_Pair.index(orig_list_key), (edge_cost_orig + delta_score))

##-----------------------------------------------------
''' this function modules the UpdateEdgeCost_Conn_Reln function
arguments: nodeA --- whose neighborhood (x) needs to be reviewed
	   nodeB --- relationship between this node and neighborhood (x) of nodeA needs to be examined
	   neighborhood_edge_type - type of neighborhood of nodeA that is to be examined
	   target_edge_type: species the cost of this type of edge between nodeB and x that needs to be updated
	   node_A_B_edge_type: edge type between nodeA and nodeB
'''
def ModuleUpdCost(Reachability_Graph_Mat, nodeA, nodeB, neighborhood_edge_type, target_edge_type, node_A_B_edge_type):    
  
  global tree_count
  
  # here an edge type between nodeA and nodeB is to be evaluated
  # define the neighborhood of node A
  nodeA_neighb = []
  if (neighborhood_edge_type == DIRECTED_IN_EDGE):
    nodeA_neighb.extend(Taxa_Info_Dict[nodeA]._GetNonProcessedInEdgeList())
  elif (neighborhood_edge_type == BI_DIRECTED_EDGE):
    nodeA_neighb.extend(Taxa_Info_Dict[nodeA]._GetNonProcessedEqEdgeList())
  elif (neighborhood_edge_type == DIRECTED_OUT_EDGE):
    nodeA_neighb.extend(Taxa_Info_Dict[nodeA]._GetNonProcessedOutEdgeList())
  else:
    nodeA_neighb.extend(Taxa_Info_Dict[nodeA]._GetNonProcessedNoEdgeList())

  if nodeB in nodeA_neighb:
    nodeA_neighb.remove(nodeB)	# dont consider nodeB    
   
  # add - sourya
  """ this section of the code checks whether any taxa is already related with nodeA,
  though the information is not reflected in the Non processed lists """
  nodeA_clust_idx = Taxa_Info_Dict[nodeA]._Get_Taxa_Part_Clust_Idx()
  nodeA_clust_reach_mat_idx = CURRENT_CLUST_IDX_LIST.index(nodeA_clust_idx)
  for x in nodeA_neighb:
    x_clust_idx = Taxa_Info_Dict[x]._Get_Taxa_Part_Clust_Idx()
    if (x_clust_idx == nodeA_clust_idx):
      nodeA_neighb.remove(x)
    else:
      x_clust_reach_mat_idx = CURRENT_CLUST_IDX_LIST.index(x_clust_idx)
      # check if there exists a connection
      if (Reachability_Graph_Mat[nodeA_clust_reach_mat_idx][x_clust_reach_mat_idx] > 0) \
	  or (Reachability_Graph_Mat[x_clust_reach_mat_idx][nodeA_clust_reach_mat_idx] > 0):
	nodeA_neighb.remove(x)
  # end add - sourya
  
  """ now process the non connected neighborhood of nodeA """
  if (len(nodeA_neighb) > 0):
    for l in nodeA_neighb:	    
      # now we calculate the score of the derived edge connection
      # (in terms of deviation from the best possible edge between these 2 nodes) 
      key1 = (nodeB, l)
      key2 = (l, nodeB)
      if (key1 not in TaxaPair_Reln_Dict) and (key2 not in TaxaPair_Reln_Dict):
	if (target_edge_type == NO_EDGE):
	  delta_score = ORIG_NO_EDGE_REMAIN_NO_EDGE
	else:
	  delta_score = ORIG_DIFF_SRC_TREE_CONNECT_SCORE
      else:
	if key1 in TaxaPair_Reln_Dict:
	  delta_score = TaxaPair_Reln_Dict[key1]._GetConnPrVal(target_edge_type)
	  if (target_edge_type != NO_EDGE):
	    # check whether originally between these 2 nodes, NO_EDGE connection was predominant 
	    if (TaxaPair_Reln_Dict[key1]._GetConnPrVal(NO_EDGE) > 0):
	      delta_score = delta_score + ORIG_NO_EDGE_BECOME_CONNECTED
	else:
	  # at first complement the target edge type
	  if (target_edge_type == DIRECTED_IN_EDGE):
	    delta_score = TaxaPair_Reln_Dict[key2]._GetConnPrVal(DIRECTED_OUT_EDGE)
	  elif (target_edge_type == DIRECTED_OUT_EDGE):
	    delta_score = TaxaPair_Reln_Dict[key2]._GetConnPrVal(DIRECTED_IN_EDGE)
	  else:
	    delta_score = TaxaPair_Reln_Dict[key2]._GetConnPrVal(target_edge_type)
	  if (target_edge_type != NO_EDGE):
	    # check whether originally between these 2 nodes, NO_EDGE connection was predominant 
	    if (TaxaPair_Reln_Dict[key2]._GetConnPrVal(NO_EDGE) > 0):
	      delta_score = delta_score + ORIG_NO_EDGE_BECOME_CONNECTED
	      
      # here we update the score of edge between nodeA and l 
      # (where l belongs to the neighborhood of nodeA)
      # now call the function to update the delta_score of list element
      """
      Important - sourya
      adjustment of cost metrics on Cost_List_Taxa_Pair_Single_Reln is commented for the moment
      as it produces better performance 
      """
      kt1 = (nodeA, l)
      kt2 = (l, nodeA)
      if kt1 in TaxaPair_Reln_Dict:
	#UpdateEdgeCostInMemList(Cost_List_Taxa_Pair_Single_Reln, nodeA, l, delta_score, neighborhood_edge_type)
	UpdateEdgeCostInMemList(Cost_List_Taxa_Pair_Multi_Reln, nodeA, l, delta_score, neighborhood_edge_type)
      elif kt2 in TaxaPair_Reln_Dict:
	# at first complement the neighborhood_edge_type
	if (neighborhood_edge_type == DIRECTED_IN_EDGE):
	  #UpdateEdgeCostInMemList(Cost_List_Taxa_Pair_Single_Reln, l, nodeA, delta_score, DIRECTED_OUT_EDGE)  
	  UpdateEdgeCostInMemList(Cost_List_Taxa_Pair_Multi_Reln, l, nodeA, delta_score, DIRECTED_OUT_EDGE)  
	elif (neighborhood_edge_type == DIRECTED_OUT_EDGE):
	  #UpdateEdgeCostInMemList(Cost_List_Taxa_Pair_Single_Reln, l, nodeA, delta_score, DIRECTED_IN_EDGE)
	  UpdateEdgeCostInMemList(Cost_List_Taxa_Pair_Multi_Reln, l, nodeA, delta_score, DIRECTED_IN_EDGE)
	else:
	  #UpdateEdgeCostInMemList(Cost_List_Taxa_Pair_Single_Reln, l, nodeA, delta_score, neighborhood_edge_type)
	  UpdateEdgeCostInMemList(Cost_List_Taxa_Pair_Multi_Reln, l, nodeA, delta_score, neighborhood_edge_type)

##-----------------------------------------------------
''' this function is called to update the edge score between a set of taxa
nodeA and nodeB are two input nodes, which are already connected in the final tree
i.e. their final edge type is decided

here the NON-PROCESSED neighborhood of nodeA and nodeB are checked

1)(say x belongs to non processed neighborhood of nodeB, with ? as the edge type) --- 
  check the possible connection of x to nodeA
  the score of this connection from x to nodeA is the priority associated with corresponding edge type
  this score is added with the score of edge type ? between nodeB to x
2)(say y belongs to non processed neighborhood of nodeA, with ? as the edge type) --- 
  check the possible connection of y to nodeB
  the score of this connection from y to nodeB is the priority associated with corresponding edge type
  this score is added with the score of edge type ? between nodeA to y '''
  
def UpdateEdgeCost_Conn_Reln(Reachability_Graph_Mat, nodeA, nodeB, edge_type):
  
  # case 1 - A->B connection 
  if (edge_type == DIRECTED_OUT_EDGE):        
    ##-----------------------------------------
    ## check non processed neighborhood of nodeA
    ##-----------------------------------------
    # check the non processed neighbors of nodeA that can be connected in future via edge C->A
    # check the possible connection C->B 
    ModuleUpdCost(Reachability_Graph_Mat, nodeA, nodeB, DIRECTED_IN_EDGE, DIRECTED_IN_EDGE, edge_type)    
    # check the non processed neighbors of nodeA that can be connected in future via edge C=A
    # check the possible connection C->B
    ModuleUpdCost(Reachability_Graph_Mat, nodeA, nodeB, BI_DIRECTED_EDGE, DIRECTED_IN_EDGE, edge_type)
    # add - sourya
    # check the non processed neighbors of nodeA that can be connected in future via edge C><A
    # check the possible connection C><B
    ModuleUpdCost(Reachability_Graph_Mat, nodeA, nodeB, NO_EDGE, NO_EDGE, edge_type)
    # end add - sourya
    ##-----------------------------------------
    ## check non processed neighborhood of nodeB
    ##-----------------------------------------
    # check the non processed neighbors of node B that can be connected in future via edge B->C
    # check the possible connection A->C
    ModuleUpdCost(Reachability_Graph_Mat, nodeB, nodeA, DIRECTED_OUT_EDGE, DIRECTED_OUT_EDGE, edge_type)    
    # check the non processed neighbors of node B that can be connected in future via edge B=C
    # check the possible connection A->C 
    ModuleUpdCost(Reachability_Graph_Mat, nodeB, nodeA, BI_DIRECTED_EDGE, DIRECTED_OUT_EDGE, edge_type)    
      
  # case 2 - A<-B connection 
  if (edge_type == DIRECTED_IN_EDGE):    
    ##-----------------------------------------
    ## check non processed neighborhood of nodeA 
    ##-----------------------------------------
    # check the non processed neighbors of node A that can be connected in future via edge A->C
    # check the possible connection B->C
    ModuleUpdCost(Reachability_Graph_Mat, nodeA, nodeB, DIRECTED_OUT_EDGE, DIRECTED_OUT_EDGE, edge_type)    
    # check the non processed neighbors of node A that can be connected in future via edge A=C
    # check the possible connection B->C 
    ModuleUpdCost(Reachability_Graph_Mat, nodeA, nodeB, BI_DIRECTED_EDGE, DIRECTED_OUT_EDGE, edge_type)
    ##-----------------------------------------
    ## check non processed neighborhood of nodeB
    ##-----------------------------------------
    # check the non processed neighbors of node B that can be connected in future via edge C->B
    # establish the possible connection C->A
    ModuleUpdCost(Reachability_Graph_Mat, nodeB, nodeA, DIRECTED_IN_EDGE, DIRECTED_IN_EDGE, edge_type)
    # check the neighbors of node B that can be connected in future via edge C=B
    # establish the connection C->A
    ModuleUpdCost(Reachability_Graph_Mat, nodeB, nodeA, BI_DIRECTED_EDGE, DIRECTED_IN_EDGE, edge_type)    		    
    # add - sourya
    # check the non processed neighbors of nodeB that can be connected in future via edge C><B
    # check the possible connection C><A 
    ModuleUpdCost(Reachability_Graph_Mat, nodeB, nodeA, NO_EDGE, NO_EDGE, edge_type)
    # end add - sourya 
    
  # case 3 - A=B connection 
  if (edge_type == BI_DIRECTED_EDGE):    
    ##-----------------------------------------
    ## check non processed neighborhood of nodeA
    ##-----------------------------------------
    # check the non processed neighbors of node A that can be connected in future via edge A->C
    # check the possible connection B->C
    ModuleUpdCost(Reachability_Graph_Mat, nodeA, nodeB, DIRECTED_OUT_EDGE, DIRECTED_OUT_EDGE, edge_type)
    # check the non processed neighbors of node A that can be connected in future via edge C->A 
    # establish the connection C->B
    ModuleUpdCost(Reachability_Graph_Mat, nodeA, nodeB, DIRECTED_IN_EDGE, DIRECTED_IN_EDGE, edge_type)
    # check the non processed neighbors of node A that can be connected in future via edge A=C
    # check the possible connection B=C
    ModuleUpdCost(Reachability_Graph_Mat, nodeA, nodeB, BI_DIRECTED_EDGE, BI_DIRECTED_EDGE, edge_type)
    # add - sourya
    # check the non processed neighbors of node A that can be connected in future via edge A><C
    # check the possible connection B><C
    ModuleUpdCost(Reachability_Graph_Mat, nodeA, nodeB, NO_EDGE, NO_EDGE, edge_type)
    # end add - sourya
    ##-----------------------------------------
    ## check non processed neighborhood of nodeB 
    ##-----------------------------------------
    # check the non processed neighbors of node B that can be connected in future via edge B->C
    # check the connection A->C
    ModuleUpdCost(Reachability_Graph_Mat, nodeB, nodeA, DIRECTED_OUT_EDGE, DIRECTED_OUT_EDGE, edge_type)
    # check the non processed neighbors of node B that can be connected in future via edge C->B 
    # establish the possible connection C->A
    ModuleUpdCost(Reachability_Graph_Mat, nodeB, nodeA, DIRECTED_IN_EDGE, DIRECTED_IN_EDGE, edge_type)    
    # check the non processed neighbors of node B that can be connected in future via edge B=C
    # check the possible connection A=C
    ModuleUpdCost(Reachability_Graph_Mat, nodeB, nodeA, BI_DIRECTED_EDGE, BI_DIRECTED_EDGE, edge_type)    
    # add - sourya
    # check the non processed neighbors of nodeB that can be connected in future via edge B><C
    # check the possible connection A><C  
    ModuleUpdCost(Reachability_Graph_Mat, nodeB, nodeA, NO_EDGE, NO_EDGE, edge_type)
    # end add - sourya
		      
  # case 4 - A><B connection (no connection)
  if (edge_type == NO_EDGE):    
    ##-----------------------------------------
    ## check non processed neighborhood of nodeA 
    ##-----------------------------------------
    # check the non processed neighbors of node A that can be connected in future via edge A->C
    # check the possible connection B><C
    ModuleUpdCost(Reachability_Graph_Mat, nodeA, nodeB, DIRECTED_OUT_EDGE, NO_EDGE, edge_type)
    # check the non processed neighbors of node A that can be connected in future via edge A=C
    # check the possible connection B><C      
    ModuleUpdCost(Reachability_Graph_Mat, nodeA, nodeB, BI_DIRECTED_EDGE, NO_EDGE, edge_type)
    ##-----------------------------------------
    ## check non processed neighborhood of nodeB 
    ##-----------------------------------------
    # check the non processed neighbors of node B that can be connected in future via edge B->C
    # check the possible connection A><C
    ModuleUpdCost(Reachability_Graph_Mat, nodeB, nodeA, DIRECTED_OUT_EDGE, NO_EDGE, edge_type)
    # check the non processed neighbors of node B that can be connected in future via edge B=C
    # check the possible connection A><C
    ModuleUpdCost(Reachability_Graph_Mat, nodeB, nodeA, BI_DIRECTED_EDGE, NO_EDGE, edge_type)    

#------------------------------------------------------
""" this function sorts the input list containing the costs of different relations between individual taxa pairs
we don't use python standard sort routine
rather we use sorting giving different weightage to different edge types, in case of identical cost values """
def Sort_Cost_List_Initial(Cost_List_Node_Pair):
  Sort_Priority_Queue(Cost_List_Node_Pair)
  
#------------------------------------------------------
""" this code section implements the max priority queue """
#------------------------------------------------------
# parent node of current node
def Parent(idx):
  return int((idx-1) / 2)

# left child node of current node  
def Left(idx):
  return (2*idx+1)

# right child node of current node
def Right(idx):
  return (2*idx+2)

#------------------------------------------------------------------------
# version 4 - latest - used for the final implementation
#-------------------------
"""
def Lower_Score_Value(inp_queue, i, j):
  #---------------------
  # two cases 
  #---------------------
  # case 1 - when both of the scores are non-positive
  if (inp_queue[i][3] <= 0) and (inp_queue[j][3] <= 0):
    if (TaxaPair_Reln_Dict[(inp_queue[i][0], inp_queue[i][1])]._GetEdgeWeight(inp_queue[i][2])\
	  < TaxaPair_Reln_Dict[(inp_queue[j][0], inp_queue[j][1])]._GetEdgeWeight(inp_queue[j][2])):
      return 1      
    elif (TaxaPair_Reln_Dict[(inp_queue[j][0], inp_queue[j][1])]._GetEdgeWeight(inp_queue[j][2])\
	  < TaxaPair_Reln_Dict[(inp_queue[i][0], inp_queue[i][1])]._GetEdgeWeight(inp_queue[i][2])):
      return 0      
    elif (inp_queue[i][3] < inp_queue[j][3]):
      return 1
    elif (inp_queue[i][3] > inp_queue[j][3]):
      return 0
    elif ((inp_queue[i][2] == NO_EDGE) and (inp_queue[j][2] != NO_EDGE)):
      # no edge has high priority
      return 0
    elif ((inp_queue[j][2] == NO_EDGE) and (inp_queue[i][2] != NO_EDGE)):
      return 1
    elif ((inp_queue[i][2] == BI_DIRECTED_EDGE) and ((inp_queue[j][2] == DIRECTED_OUT_EDGE) \
      or (inp_queue[j][2] == DIRECTED_IN_EDGE))):
      return 1      
    elif ((inp_queue[j][2] == BI_DIRECTED_EDGE) and ((inp_queue[i][2] == DIRECTED_OUT_EDGE) \
      or (inp_queue[i][2] == DIRECTED_IN_EDGE))):
      return 0
    elif (TaxaPair_Reln_Dict[(inp_queue[i][0], inp_queue[i][1])]._GetConnPrVal(inp_queue[i][2])\
	  < TaxaPair_Reln_Dict[(inp_queue[j][0], inp_queue[j][1])]._GetConnPrVal(inp_queue[j][2])):
      return 1
    elif (TaxaPair_Reln_Dict[(inp_queue[j][0], inp_queue[j][1])]._GetConnPrVal(inp_queue[j][2])\
	  < TaxaPair_Reln_Dict[(inp_queue[i][0], inp_queue[i][1])]._GetConnPrVal(inp_queue[i][2])):
      return 0    
  else:
    # case 2 - otherwise we follow previous algorithm 
    if (inp_queue[i][3] < inp_queue[j][3]):
      return 1
    elif (inp_queue[i][3] == inp_queue[j][3]):
      # for tie case of cost
      if ((inp_queue[i][2] == NO_EDGE) and (inp_queue[j][2] != NO_EDGE)):
	# change - sourya
	# no edge has high priority
	return 0
      elif ((inp_queue[j][2] == NO_EDGE) and (inp_queue[i][2] != NO_EDGE)):
	return 1
      # end change - sourya
      elif ((inp_queue[i][2] == BI_DIRECTED_EDGE) and ((inp_queue[j][2] == DIRECTED_OUT_EDGE) or (inp_queue[j][2] == DIRECTED_IN_EDGE))):
	# bi directed edge is set to low priority, compared to other directed edges - add - sourya
	return 1      
      elif ((inp_queue[j][2] == BI_DIRECTED_EDGE) and ((inp_queue[i][2] == DIRECTED_OUT_EDGE) or (inp_queue[i][2] == DIRECTED_IN_EDGE))):
	return 0
      elif (TaxaPair_Reln_Dict[(inp_queue[i][0], inp_queue[i][1])]._GetConnPrVal(inp_queue[i][2])\
	    < TaxaPair_Reln_Dict[(inp_queue[j][0], inp_queue[j][1])]._GetConnPrVal(inp_queue[j][2])):
	# higher edge priority have greater priority
	return 1
      elif (TaxaPair_Reln_Dict[(inp_queue[j][0], inp_queue[j][1])]._GetConnPrVal(inp_queue[j][2])\
	    < TaxaPair_Reln_Dict[(inp_queue[i][0], inp_queue[i][1])]._GetConnPrVal(inp_queue[i][2])):
	return 0    
      elif (TaxaPair_Reln_Dict[(inp_queue[i][0], inp_queue[i][1])]._GetEdgeWeight(inp_queue[i][2])\
	    < TaxaPair_Reln_Dict[(inp_queue[j][0], inp_queue[j][1])]._GetEdgeWeight(inp_queue[j][2])):
	# higher edge frequency have greater priority - add - sourya
	return 1      
      elif (TaxaPair_Reln_Dict[(inp_queue[j][0], inp_queue[j][1])]._GetEdgeWeight(inp_queue[j][2])\
	    < TaxaPair_Reln_Dict[(inp_queue[i][0], inp_queue[i][1])]._GetEdgeWeight(inp_queue[i][2])):
	return 0      
      
  return 0
"""
#------------------------------------------------------------------------
# version 3
#-------------------------
"""
def Lower_Score_Value(inp_queue, i, j):
  if (inp_queue[i][3] < inp_queue[j][3]):
    return 1
  elif (inp_queue[i][3] == inp_queue[j][3]):
    # for tie case of cost
    if ((inp_queue[i][2] == NO_EDGE) and (inp_queue[j][2] != NO_EDGE)):
      # change - sourya
      # no edge has high priority
      return 0
    elif ((inp_queue[j][2] == NO_EDGE) and (inp_queue[i][2] != NO_EDGE)):	#added
      return 1
    # end change - sourya
    elif ((inp_queue[i][2] == BI_DIRECTED_EDGE) and ((inp_queue[j][2] == DIRECTED_OUT_EDGE) or (inp_queue[j][2] == DIRECTED_IN_EDGE))):
      # bi directed edge is set to low priority, compared to other directed edges - add - sourya
      return 1      
    elif ((inp_queue[j][2] == BI_DIRECTED_EDGE) and ((inp_queue[i][2] == DIRECTED_OUT_EDGE) or (inp_queue[i][2] == DIRECTED_IN_EDGE))):	#added
      return 0
    elif (TaxaPair_Reln_Dict[(inp_queue[i][0], inp_queue[i][1])]._GetConnPrVal(inp_queue[i][2])\
	  < TaxaPair_Reln_Dict[(inp_queue[j][0], inp_queue[j][1])]._GetConnPrVal(inp_queue[j][2])):
      # higher edge priority have greater priority
      return 1
    elif (TaxaPair_Reln_Dict[(inp_queue[j][0], inp_queue[j][1])]._GetConnPrVal(inp_queue[j][2])\
	  < TaxaPair_Reln_Dict[(inp_queue[i][0], inp_queue[i][1])]._GetConnPrVal(inp_queue[i][2])):	#added
      return 0    
    elif (TaxaPair_Reln_Dict[(inp_queue[i][0], inp_queue[i][1])]._GetEdgeWeight(inp_queue[i][2])\
	  < TaxaPair_Reln_Dict[(inp_queue[j][0], inp_queue[j][1])]._GetEdgeWeight(inp_queue[j][2])):
      # higher edge frequency have greater priority - add - sourya
      return 1      
    elif (TaxaPair_Reln_Dict[(inp_queue[j][0], inp_queue[j][1])]._GetEdgeWeight(inp_queue[j][2])\
	  < TaxaPair_Reln_Dict[(inp_queue[i][0], inp_queue[i][1])]._GetEdgeWeight(inp_queue[i][2])):	#added
      return 0      
    
  return 0
"""
#------------------------------------------------------------------------
# version 2
#-------------------------
# this function is used for cost based sorting of the inp_queue in the unweighted supertree
# this function compares two elements of the heap
# returns 1 if element in the i'th index is lower than element in the j'th index
# that is, j'th index has higher priority than i'th index
"""
def Lower_Score_Value(inp_queue, i, j):
  if (inp_queue[i][3] < inp_queue[j][3]):
    return 1
  elif (inp_queue[i][3] == inp_queue[j][3]):
    # following checks are performed in tie case of cost
    if (((inp_queue[i][2] == DIRECTED_OUT_EDGE) or (inp_queue[i][2] == DIRECTED_IN_EDGE)) and \
	((inp_queue[j][2] == BI_DIRECTED_EDGE) or (inp_queue[j][2] == NO_EDGE))):
      return 0
    elif (((inp_queue[j][2] == DIRECTED_OUT_EDGE) or (inp_queue[j][2] == DIRECTED_IN_EDGE)) and \
	((inp_queue[i][2] == BI_DIRECTED_EDGE) or (inp_queue[i][2] == NO_EDGE))):
      return 1
    elif ((inp_queue[i][2] == NO_EDGE) and (inp_queue[j][2] == BI_DIRECTED_EDGE)):
      return 0
    elif ((inp_queue[j][2] == NO_EDGE) and (inp_queue[i][2] == BI_DIRECTED_EDGE)):
      return 1
    elif (TaxaPair_Reln_Dict[(inp_queue[i][0], inp_queue[i][1])]._GetConnPrVal(inp_queue[i][2])\
	  < TaxaPair_Reln_Dict[(inp_queue[j][0], inp_queue[j][1])]._GetConnPrVal(inp_queue[j][2])):
      # higher edge priority have greater priority
      return 1
    elif (TaxaPair_Reln_Dict[(inp_queue[i][0], inp_queue[i][1])]._GetConnPrVal(inp_queue[i][2])\
	  > TaxaPair_Reln_Dict[(inp_queue[j][0], inp_queue[j][1])]._GetConnPrVal(inp_queue[j][2])):
      # higher edge priority have greater priority
      return 0    
    elif (TaxaPair_Reln_Dict[(inp_queue[i][0], inp_queue[i][1])]._GetEdgeWeight(inp_queue[i][2])\
	  < TaxaPair_Reln_Dict[(inp_queue[j][0], inp_queue[j][1])]._GetEdgeWeight(inp_queue[j][2])):
      # higher edge frequency have greater priority - add - sourya
      return 1      
    elif (TaxaPair_Reln_Dict[(inp_queue[i][0], inp_queue[i][1])]._GetEdgeWeight(inp_queue[i][2])\
	  > TaxaPair_Reln_Dict[(inp_queue[j][0], inp_queue[j][1])]._GetEdgeWeight(inp_queue[j][2])):
      # higher edge frequency have greater priority - add - sourya
      return 0          
  return 0
"""
#------------------------------------------------------------------------
# version 1 
#-------------------------

# this function is used for cost based sorting of the inp_queue in the unweighted supertree
# this function compares two elements of the heap
# returns 1 if element in the i'th index is lower than element in the j'th index

def Lower_Score_Value(inp_queue, i, j):
  if (inp_queue[i][3] < inp_queue[j][3]):
    return 1
  elif (inp_queue[i][3] > inp_queue[j][3]):
    return 0
  else:	#if (inp_queue[i][3] == inp_queue[j][3]):
    # for tie case of cost
    if ((inp_queue[i][2] == NO_EDGE) and (inp_queue[j][2] != NO_EDGE)):
      # no edge has low priority
      return 1
    elif ((inp_queue[j][2] == NO_EDGE) and (inp_queue[i][2] != NO_EDGE)):	#added
      return 0
    elif ((inp_queue[i][2] == BI_DIRECTED_EDGE) and ((inp_queue[j][2] == DIRECTED_OUT_EDGE) or (inp_queue[j][2] == DIRECTED_IN_EDGE))):
      # bi directed edge is set to low priority, compared to other directed edges - add - sourya
      return 1      
    elif ((inp_queue[j][2] == BI_DIRECTED_EDGE) and ((inp_queue[i][2] == DIRECTED_OUT_EDGE) or (inp_queue[i][2] == DIRECTED_IN_EDGE))):	#added
      return 0
    elif (TaxaPair_Reln_Dict[(inp_queue[i][0], inp_queue[i][1])]._GetConnPrVal(inp_queue[i][2])\
	  < TaxaPair_Reln_Dict[(inp_queue[j][0], inp_queue[j][1])]._GetConnPrVal(inp_queue[j][2])):
      # higher edge priority have greater priority
      return 1
    elif (TaxaPair_Reln_Dict[(inp_queue[j][0], inp_queue[j][1])]._GetConnPrVal(inp_queue[j][2])\
	  < TaxaPair_Reln_Dict[(inp_queue[i][0], inp_queue[i][1])]._GetConnPrVal(inp_queue[i][2])):	#added
      return 0    
    elif (TaxaPair_Reln_Dict[(inp_queue[i][0], inp_queue[i][1])]._GetEdgeWeight(inp_queue[i][2])\
	  < TaxaPair_Reln_Dict[(inp_queue[j][0], inp_queue[j][1])]._GetEdgeWeight(inp_queue[j][2])):
      # higher edge frequency have greater priority - add - sourya
      return 1      
    elif (TaxaPair_Reln_Dict[(inp_queue[j][0], inp_queue[j][1])]._GetEdgeWeight(inp_queue[j][2])\
	  < TaxaPair_Reln_Dict[(inp_queue[i][0], inp_queue[i][1])]._GetEdgeWeight(inp_queue[i][2])):	#added
      return 0      
    
  return 0

#------------------------------------------------------------------------
  
# this function exchanges two elements in the heap
def Exchange_Elem(inp_queue, i, j):
  temp_key = inp_queue[i][0]
  inp_queue[i][0] = inp_queue[j][0]
  inp_queue[j][0] = temp_key
  temp_key = inp_queue[i][1]
  inp_queue[i][1] = inp_queue[j][1]
  inp_queue[j][1] = temp_key
  temp_edge = inp_queue[i][2]
  inp_queue[i][2] = inp_queue[j][2]
  inp_queue[j][2] = temp_edge
  temp_val = inp_queue[i][3]
  inp_queue[i][3] = inp_queue[j][3]
  inp_queue[j][3] = temp_val

# maintain max heap property
# note: heap_size may not be the actual length of the queue
# but the working length (on which the remaining sorting operation will take place)
def Max_Heapify(inp_queue, idx, heap_size):
  l = Left(idx)
  r = Right(idx)
  if (l < heap_size) and Lower_Score_Value(inp_queue, idx, l):
    largest_idx = l
  else:
    largest_idx = idx
  if (r < heap_size) and Lower_Score_Value(inp_queue, largest_idx, r):
    largest_idx = r
  if (largest_idx != idx):
    Exchange_Elem(inp_queue, idx, largest_idx)
    Max_Heapify(inp_queue, largest_idx, heap_size)

# extract the current maximum and also pop the element from the heap
def Heap_Extract_Max(inp_queue):
  if (len(inp_queue) < 1):
    print 'underflow of max priority queue'
  # 1st element is the maximum
  max_elem = list(inp_queue[0])
  # replace the first element with the last element of the queue
  Exchange_Elem(inp_queue, 0, len(inp_queue) - 1)
  # delete the last element of the queue
  del inp_queue[len(inp_queue) - 1]
  heap_size = len(inp_queue)
  # call the max_heapify function to maintain the heap property
  # 0 is the starting index of the list storing the heap structure
  Max_Heapify(inp_queue, 0, heap_size)
  return max_elem
  
# this function builds the priority queue (max heap property)
def Build_Max_Heap(inp_queue):
  heap_size = len(inp_queue)
  for idx in range(int(len(inp_queue) / 2), -1, -1):
    Max_Heapify(inp_queue, idx, heap_size)
    
# this is the heap sort algorithm
def Sort_Priority_Queue(inp_queue):
  Build_Max_Heap(inp_queue)
  heap_size = len(inp_queue)
  # this code is not required - sourya
  # all the sorting related stuff iteratively has been performed in the function proc_queue
  #for idx in range((len(inp_queue) - 1), 0, -1):
    #Exchange_Elem(inp_queue, 0, idx)
    #heap_size = heap_size - 1
    #Max_Heapify(inp_queue, idx, heap_size)
  # end comment - sourya
  
def Queue_Increase_Score(inp_queue, idx, target_val):
  inp_queue[idx][3] = target_val
  while (idx > 0) and Lower_Score_Value(inp_queue, Parent(idx), idx):
    Exchange_Elem(inp_queue, idx, Parent(idx))
    idx = Parent(idx)
    
def Queue_Decrease_Score(inp_queue, idx, target_val):
  inp_queue[idx][3] = target_val
  Max_Heapify(inp_queue, idx, len(inp_queue))
      