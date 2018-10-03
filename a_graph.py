#!/usr/bin/env python

######################################################################################
##	a_graph.py
#
#	Author: Michele Berselli
#		University of Padova
#		berselli.michele@gmail.com
#
##	LICENSE:
#		Copyright (C) 2018 Michele Berselli
#
#		This program is free software: you can redistribute it and/or modify
#		it under the terms of the GNU General Public License as published by
#		the Free Software Foundation.
#
#		This program is distributed in the hope that it will be useful,
#		but WITHOUT ANY WARRANTY; without even the implied warranty of
#		MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#		GNU General Public License for more details.
#
#		You should have received a copy of the GNU General Public License
#		along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
######################################################################################


########################################
####		CLASS graph             ####
########################################

class graph(object):
	''' a class to implement a basic graph '''

	class node(object):
		''' a class to implement a basic node '''

		def __init__(self, name):
			self.__name = name	#name of the node
			self.__nodi = []	#linked nodes
		#end def __init__

		def get_name(self):
			return self.__name
		#end def get_name

		def add_node(self, node_to_add):
			self.__nodi.append(node_to_add)
		#end def add_node

		def get_nodi(self):
			return self.__nodi
		#end def get_nodi

	#end class node

	def __init__(self, root):
		self.__root = graph.node(root)
	#end def __init__

	def get_root(self):
		return self.__root
	#end def get_root

	def add_node(self, node, node_to_add):
		node.add_node(node_to_add)
	#end def add_node

	########################################
	####		DFS exploration         ####
	########################################

	def explore_DFS(self, node):
		''' explore a graph in depth first mode,
		explore the graph from the given node (recursive) '''
		results = []
		self.__routine_explore_DFS(node, results)
		return results
	#end def explore_DFS

	def __routine_explore_DFS(self, node, results):
		if not node:
			return
		#end if
		[self.__routine_explore_DFS(node_i, results) for node_i in node.get_nodi()]
		results.append(node)
	#end def __routine_explore_DFS

	########################################
	####		DFS print paths         ####
	########################################

	def get_paths(self, node, length):
		''' return all the paths long length
		starting from the given node (recursive) '''
		out, results, current = [], [], 0
		self.__routine_get_paths(node, current, length, out, results)
		return out
	#end def get_paths

	def __routine_get_paths(self, node, current, length, out, results):
		if not node:
			return
		else:
			current += 1
			results.append(node)
		#end if
		if not length == current:
			[self.__routine_get_paths(node_i, current, length, out, results) for node_i in node.get_nodi()]
		else:
			out.append(results[:])
		#end if
		results.pop()
	#end def __routine_get_paths

	########################################
	####		Longest Path            ####
	######################################## 

	def longest_path(self, node):
		''' explore a graph in depth first mode,
		explore the graph from the given node and calculate 
		longest path to leafs from the starting node (recursive) '''
		longest, current = [0], 0
		self.__routine_longest_path(node, current, longest)
		return longest[0]
	#end def longest_path

	def __routine_longest_path(self, node, current, longest):
		if not node:
			return
		else:
			current += 1
			if longest[0] < current:
				longest[0] = current
			else:
				return
			#end if
		#end if
		[self.__routine_longest_path(node_i, current, longest) for node_i in node.get_nodi()]
	#end def __routine_longest_path

	########################################
	####		Path Check              ####
	########################################

	def path_check(self, node, length):
		''' explore a graph in depth first mode,
		explore the tree from the given node and check if
		there is a path at least long length (recursive) '''
		current, check = 0, [0]
		self.__routine_path_check(node, current, length, check)
		return bool(check[0])
	#end def path_check

	def __routine_path_check(self, node, current, length, check):
		if not node:
			current -= 1
			return
		else:
			current += 1
			if check[0]:
				return
			elif length == current:
				check[0] = 1
				return
			#end if
		#end if
		[self.__routine_path_check(node_i, current, length, check) for node_i in node.get_nodi()]
	#end def __routine_path_check

#end graph


########################################
####		TESTING                 ####
########################################

#a_graph = graph(0)
#a_node = graph.node(1)
#b_node = graph.node(2)
#c_node = graph.node(3)
#d_node = graph.node(4)
#e_node = graph.node(5)
#f_node = graph.node(6)
#h_node = graph.node(7)

#a_graph.add_node(a_node, b_node)
#a_graph.add_node(a_node, c_node)
#a_graph.add_node(b_node, c_node)
#a_graph.add_node(b_node, d_node)
#a_graph.add_node(b_node, e_node)
#a_graph.add_node(c_node, e_node)
#a_graph.add_node(c_node, f_node)
#a_graph.add_node(d_node, f_node)
#a_graph.add_node(e_node, f_node)
#a_graph.add_node(f_node, h_node)

#	      6 - 7 
#	     /
#	    3 - 5 - 6 - 7
#	   /
#	  2 - 4 - 6 - 7
#	 / \
#	1   5 - 6 - 7
#	 \
#	  3 - 5 - 6 - 7
#	   \ 
#	    6 - 7





