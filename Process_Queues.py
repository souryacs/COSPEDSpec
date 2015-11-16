#!/usr/bin/env python

import Header
from Header import *
import Cost_Update
from Cost_Update import *
import ReachGraph_Update
from ReachGraph_Update import *
import Conflict_Detect
from Conflict_Detect import *

#-------------------------------------------------------
""" 
this function processes individual queues (score lists) to construct the final supertree 

variables: 
Single_Reln_no_Conflict_Queue_process: 
if 1, then points to the non-conflicting queue 
otherwise, points to the multi relation (conflicting) queue
"""
def Proc_Queue(Reachability_Graph_Mat, Single_Reln_no_Conflict_Queue_process, Output_Text_File):
	if (Single_Reln_no_Conflict_Queue_process == 1):
		Inp_Queue = Cost_List_Taxa_Pair_Single_Reln
	else:
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
			if (Single_Reln_no_Conflict_Queue_process == 1):
				fp.write('\n ===>> NON CONFLICTING QUEUE -- ')
			else:
				fp.write('\n ===>> CONFLICTING QUEUE -- ')      
			fp.write(' current extracted max element: ' + str(src_taxa_label) + ' and ' + str(dest_taxa_label) + \
								' reln type: ' + str(reln_type) + ' conn score: ' + str(conn_score))
			fp.close()
		
		""" 
		if the current support score based relation does not induce a conflict to the existing configuration of the final supertree
		then include the current relation (and resolve corresponding couplet) in it 
		"""
		conflict_detection = Possible_Conflict_Curr_Reln(src_taxa_label, dest_taxa_label, Reachability_Graph_Mat, reln_type, Output_Text_File)
		
		if (conflict_detection == 0):
			""" 
			current element does not create a cycle / conflict
			not that it is already present in the supertree 
			valid connection is found - append this connection to the final formed tree 
			"""
			if (DEBUG_LEVEL >= 2):
				fp = open(Output_Text_File, 'a')    
				if (Single_Reln_no_Conflict_Queue_process == 1):
					queue_str = 'NON CONFLICTING QUEUE (higher priority)'
				else:
					queue_str = 'CONFLICTING QUEUE (lower priority)'
				fp.write('\n ==>>>>>>>>> NEW CONN --- ' + str(queue_str) + 'nodes to be connected: ' + str(src_taxa_label) + ' and ' + str(dest_taxa_label) + \
							' nodes indices: ' + str(COMPLETE_INPUT_TAXA_LIST.index(src_taxa_label)) + ' and ' + str(COMPLETE_INPUT_TAXA_LIST.index(dest_taxa_label)) + \
							' edge type: ' + str(reln_type) + ' conn score: ' + str(conn_score))
				fp.close()
			
			# also update the reachability graph information
			Reachability_Graph_Mat = AdjustReachGraph(Reachability_Graph_Mat, src_taxa_label, dest_taxa_label, reln_type)
									
	return Reachability_Graph_Mat
  