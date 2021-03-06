#!/usr/bin/env python

import Header
from Header import *
  
#------------------------------------------------------
""" 
this code section implements the max priority queue 
"""
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
# version 2 - this is the latest - sourya - 1st nov 2015
#-------------------------
"""
this function is used for cost based sorting of the inp_queue in the unweighted supertree
this function compares two elements of the heap
returns 1 if element in the i'th index is lower than element in the j'th index
that is, j'th index has higher priority than i'th index
"""
def Lower_Score_Value(inp_queue, i, j):
	key1 = (inp_queue[i][0], inp_queue[i][1])
	reln1 = inp_queue[i][2]
	score1 = inp_queue[i][3]

	key2 = (inp_queue[j][0], inp_queue[j][1])
	reln2 = inp_queue[j][2]
	score2 = inp_queue[j][3]

	if (score1 > 0) or (score2 > 0):
		if (score1 < score2):
			return 1
		elif (score1 > score2):
			return 0
		elif (score1 == score2):
			# following checks are performed in tie case of cost
			if (((reln1 == RELATION_R1) or (reln1 == RELATION_R2)) and ((reln2 == RELATION_R3) or (reln2 == RELATION_R4))):
				return 0
			elif (((reln2 == RELATION_R1) or (reln2 == RELATION_R2)) and ((reln1 == RELATION_R3) or (reln1 == RELATION_R4))):
				return 1
			elif ((reln1 == RELATION_R4) and (reln2 == RELATION_R3)):
				return 0
			elif ((reln2 == RELATION_R4) and (reln1 == RELATION_R3)):
				return 1
			elif (TaxaPair_Reln_Dict[key1]._GetConnPrVal(reln1) < TaxaPair_Reln_Dict[key2]._GetConnPrVal(reln2)):
				return 1
			elif (TaxaPair_Reln_Dict[key1]._GetConnPrVal(reln1) > TaxaPair_Reln_Dict[key2]._GetConnPrVal(reln2)):
				return 0    
			elif (TaxaPair_Reln_Dict[key1]._GetEdgeWeight(reln1) < TaxaPair_Reln_Dict[key2]._GetEdgeWeight(reln2)):
				return 1      
			elif (TaxaPair_Reln_Dict[key1]._GetEdgeWeight(reln1) > TaxaPair_Reln_Dict[key2]._GetEdgeWeight(reln2)):
				return 0
	else:
		# here both scores are either zero or negative
		if (TaxaPair_Reln_Dict[key1]._GetEdgeWeight(reln1) < TaxaPair_Reln_Dict[key2]._GetEdgeWeight(reln2)):
			return 1      
		elif (TaxaPair_Reln_Dict[key1]._GetEdgeWeight(reln1) > TaxaPair_Reln_Dict[key2]._GetEdgeWeight(reln2)):
			return 0     
		elif (TaxaPair_Reln_Dict[key1]._GetConnPrVal(reln1) < TaxaPair_Reln_Dict[key2]._GetConnPrVal(reln2)):
			return 1
		elif (TaxaPair_Reln_Dict[key1]._GetConnPrVal(reln1) > TaxaPair_Reln_Dict[key2]._GetConnPrVal(reln2)):
			return 0    
		elif (((reln1 == RELATION_R1) or (reln1 == RELATION_R2)) and ((reln2 == RELATION_R3) or (reln2 == RELATION_R4))):
			return 0
		elif (((reln2 == RELATION_R1) or (reln2 == RELATION_R2)) and ((reln1 == RELATION_R3) or (reln1 == RELATION_R4))):
			return 1
		elif ((reln1 == RELATION_R4) and (reln2 == RELATION_R3)):
			return 0
		elif ((reln2 == RELATION_R4) and (reln1 == RELATION_R3)):
			return 1

	return 0

#--------------------------------------------------------------------
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
  
      