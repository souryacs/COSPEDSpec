#!/usr/bin/env python

import Header
from Header import *
import Cost_Update
from Cost_Update import *
import ReachGraph_Update
from ReachGraph_Update import *
import Conflict_Detect
from Conflict_Detect import *
import UtilFunc
from UtilFunc import *

#-------------------------------------------------------
""" this function processes individual queues (score lists)
to construct the final supertree 
Single_Reln_no_Conflict_Queue_process: if 1, then points to the single relation queue of score metrics
	      otherwise, points to the multi relation queue of score metrics
it can be collection of taxa pairs exhibiting single relation instance
or can be taxa pairs exhibiting multi relation instance """
def Proc_Queue(Reachability_Graph_Mat, COST_UPDATE_LATEST, Single_Reln_no_Conflict_Queue_process, \
		Output_Text_File, Bin_SupTr_Opt, Dfs_parent_refine):
  if (Single_Reln_no_Conflict_Queue_process == 1):
    Inp_Queue = Cost_List_Taxa_Pair_Single_Reln
  else:
    Inp_Queue = Cost_List_Taxa_Pair_Multi_Reln
  
  while (0 < len(Inp_Queue)):
    """ extract the 1st element of "Inp_Queue" 
    since it is sorted to have max cost at the beginning """
    outlist = Heap_Extract_Max(Inp_Queue)
    
    src_taxa_label = outlist[0]
    dest_taxa_label = outlist[1]
    src_to_dest_edge_type = outlist[2]
    conn_score = TaxaPair_Reln_Dict[(src_taxa_label, dest_taxa_label)]._GetEdgeCost_ConnReln(src_to_dest_edge_type)
    
    """ if the current score metric based relation does not induce a cycle to the existing configuration of the final supertree
    then include the connection in it """
    conflict_detection = Possible_Conflict_Curr_Reln(src_taxa_label, dest_taxa_label, Reachability_Graph_Mat, \
						      src_to_dest_edge_type, Bin_SupTr_Opt)
    
    if (conflict_detection > 0):   
      """ current element does not create a cycle / conflict
      not that it is already present in the supertree 
      valid connection is found - append this connection to the final formed tree """
      if (DEBUG_LEVEL > 0):
	fp = open(Output_Text_File, 'a')    
	if (Single_Reln_no_Conflict_Queue_process == 1):
	  queue_str = 'NON CONFLICTING QUEUE (higher priority)'
	else:
	  queue_str = 'CONFLICTING QUEUE (lower priority)'
	  
	if (conflict_detection == RELN_ALREADY_ESTABLISHED):
	  fp.write('\n **** CONFLICT -- RELN_ALREADY_ESTABLISHED ****')
	else:
	  fp.write('\n **** CONFLICT -- RELN_NOT_PERMITTED ****')
	fp.write(' ' + str(queue_str) + 'nodes involved: ' + str(src_taxa_label) + ' and ' + str(dest_taxa_label) + \
	      ' nodes indices: ' + str(COMPLETE_INPUT_TAXA_LIST.index(src_taxa_label)) + ' and ' + str(COMPLETE_INPUT_TAXA_LIST.index(dest_taxa_label)) + \
	      ' edge type: ' + str(src_to_dest_edge_type) + ' conn score: ' + str(conn_score))
	fp.close()      
      
      if (COST_UPDATE_LATEST == True):      
	""" call the selective removal function of the particular edge type and the target node type
	from the structure for individual nodes """
	if (src_to_dest_edge_type == NO_EDGE) or (src_to_dest_edge_type == BI_DIRECTED_EDGE):
	  dest_to_src_edge_type = src_to_dest_edge_type
	elif (src_to_dest_edge_type == DIRECTED_IN_EDGE):
	  dest_to_src_edge_type = DIRECTED_OUT_EDGE
	else:
	  dest_to_src_edge_type = DIRECTED_IN_EDGE
	# dest_taxa_label and src_to_dest_edge_type is removed from the structure of src_taxa_label
	Taxa_Info_Dict[src_taxa_label]._SelectivelyRemoveEdgeOrigConn(src_to_dest_edge_type, dest_taxa_label)
	# src_taxa_label and dest_to_src_edge_type is removed from the structure of dest_taxa_label
	Taxa_Info_Dict[dest_taxa_label]._SelectivelyRemoveEdgeOrigConn(dest_to_src_edge_type, src_taxa_label)
    else:
      """ current element does not create a cycle / conflict
      not that it is already present in the supertree 
      valid connection is found - append this connection to the final formed tree """
      if (DEBUG_LEVEL > 0):
	fp = open(Output_Text_File, 'a')    
	if (Single_Reln_no_Conflict_Queue_process == 1):
	  queue_str = 'NON CONFLICTING QUEUE (higher priority)'
	else:
	  queue_str = 'CONFLICTING QUEUE (lower priority)'
	fp.write('\n ==>>>>>>>>> NEW CONN --- ' + str(queue_str) + 'nodes to be connected: ' + str(src_taxa_label) + ' and ' + str(dest_taxa_label) + \
	      ' nodes indices: ' + str(COMPLETE_INPUT_TAXA_LIST.index(src_taxa_label)) + ' and ' + str(COMPLETE_INPUT_TAXA_LIST.index(dest_taxa_label)) + \
	      ' edge type: ' + str(src_to_dest_edge_type) + ' conn score: ' + str(conn_score))
	fp.close()
      
      # also update the reachability graph information
      Reachability_Graph_Mat = AdjustReachGraph(Reachability_Graph_Mat, src_taxa_label, dest_taxa_label, src_to_dest_edge_type, COST_UPDATE_LATEST)
            
      # now update the score values of the relationships transitively implied by this new relation
      if (COST_UPDATE_LATEST == True):
	for list_idx in range(len(EDGE_PROCESSED_LIST)):
	  nodeA_label = EDGE_PROCESSED_LIST[list_idx][0]
	  nodeB_label = EDGE_PROCESSED_LIST[list_idx][1]
	  edge_type = EDGE_PROCESSED_LIST[list_idx][2]
	  # now call the update edge cost function
	  UpdateEdgeCost_Conn_Reln(Reachability_Graph_Mat, nodeA_label, nodeB_label, edge_type)
	# now reset the EDGE_PROCESSED_LIST structure
	del EDGE_PROCESSED_LIST[:]
      
  return Reachability_Graph_Mat
  