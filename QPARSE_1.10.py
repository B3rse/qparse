#!/usr/bin/env python


######################################################################################
#
##	QPARSE.py (Quadruplex and Paired quAdRuplex SEarch)
#		QPARSE is a tool developed for the detection of DNA and RNA-quadruplexes forming patterns.
#		The tool is flexible and allows to perform an exhaustive search of all the possible
#		quadruplexes forming patterns of the desired base in the sequence.
#		Alternatively, patterns containing any number of islands can be detected.
#		It can also detect degenerate patterns containing imperfect islands,
#		and evaluate the symmetrical properties of the linking loops to detect longer loops
#		that can form hairpin structures.
#
#	Author: Michele Berselli
#		University of Padova
#		berselli.michele@gmail.com
#
##	LICENSE:
#		Copyright (C) 2019 Michele Berselli
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
####		LIBRARIES               ####
########################################

import sys, argparse
from a_graph import graph
import numpy as np


########################################
####		CLASS g_island          ####
########################################

class g_island(graph.node):
	''' a class to implement a graph.node object to representing a G Island '''

	def __init__(self, name, strt_idx, end_idx, G_num, gaps_num, gaps_len, score, perfect):
		graph.node.__init__(self, name)	#inheritance from graph.node
		self.__strt_idx = strt_idx	#starting index of the island
		self.__end_idx = end_idx	#ending index of the island
		self.__G_num = G_num	#number of G in the island
		self.__gaps_num = gaps_num	#number of gaps in the island
		self.__gaps_len = gaps_len	#total length of gaps in the island
		self.__score = score
		self.__perfect = perfect #bool value
	#end def __init__

	def get_strt_idx(self):
		return self.__strt_idx
	#end def get_strt_idx

	def get_end_idx(self):
		return self.__end_idx
	#end def get_end_idx

	def get_G_num(self):
		return self.__G_num
	#end def get_G_num

	def get_gaps_num(self):
		return self.__gaps_num
	#end def get_gaps_num

	def get_gaps_len(self):
		return self.__gaps_len
	#end def get_gaps_len

	def get_score(self):
		return self.__score
	#end def get_score

	def get_island(self, seq):
		return seq[self.__strt_idx:self.__end_idx + 1]
	#end def get_island

	def is_perfect(self):
		return self.__perfect
	#end def is_perfect

#end class g_Island


########################################
####		CLASS g_graph           ####
########################################

class g_graph(graph):
	''' '''

	def __init__(self, root):
		graph.__init__(self, root)	#inheritance from graph
	#end def __init__

	########################################
	####		General                 ####
	########################################

	def __routine_perfect_islands(self, quadruplex_as_list_island):
		''' return the number of perfect islands in quadruplex_as_list_island '''
		counter_perfect_islands = len(quadruplex_as_list_island)
		for island in quadruplex_as_list_island:
			if not island.is_perfect():
				counter_perfect_islands -= 1
			#end if
		#end for
		return counter_perfect_islands
	#end def __routine_perfect_islands

	def __routine_mismatch(self, quadruplex_as_list_island, min_len, bulge_only_mismatch):
		''' return if it is a mismatched PQS and the number of mismatched islands '''
		isl_dict = {} #{min_len: counts, ...}
		for island in quadruplex_as_list_island:
			isl_dict.setdefault(island.get_G_num(), 0)
			isl_dict[island.get_G_num()] += 1
		#end for
		len_isl_dict = len(isl_dict)
		if len_isl_dict == 1:
			if min_len in isl_dict:
				return False, 0
			else:
				if bulge_only_mismatch:
					is_perfect = True
					for island in quadruplex_as_list_island:
						if island.get_gaps_len() > 0:
							is_perfect = False
							break
						#end if
					#end for
					if is_perfect:
						return True, 0
					else:
						return False, 0
					#end if
				else:
					return True, 0
			#end if
		elif len_isl_dict == 2:
			return True, isl_dict[min(isl_dict)]
		else:
			return False, 0
		#end if
	#end def __routine_mismatch

	########################################
	####		DFS print paths         ####
	########################################

	def get_paths_score(self, node, length, mismatch, min_len, bulge_only_mismatch, min_perfect=1):
		''' return all the paths long length
		starting from the given node with
		an associated score '''
		out, results, current, score = [], [], 0, [0]
		self.__routine_get_paths_score(node, current, length, out, results, score, min_perfect, mismatch, min_len, bulge_only_mismatch)
		return out
	#end def get_paths_score

	def __routine_get_paths_score(self, node, current, length, out, results, score, min_perfect, mismatch, min_len, bulge_only_mismatch):
		if not node:
			return
		else:
			current += 1
			results.append(node)
			score[0] += node.get_score()
		#end if
		if not length == current:
			[self.__routine_get_paths_score(node_i, current, length, out, results, score, min_perfect, mismatch, min_len, bulge_only_mismatch) for node_i in node.get_nodi()]
		else:
			quadruplex_as_list_island = results[:]
			append_to_out = True
			if min_perfect:
				if self.__routine_perfect_islands(quadruplex_as_list_island) < min_perfect:
					append_to_out = False
				#end if
			#end if
			if mismatch:
				is_mismatch, mismatch_num = self.__routine_mismatch(quadruplex_as_list_island, min_len, bulge_only_mismatch)
				if (not is_mismatch) or (mismatch_num > mismatch):
					append_to_out = False
				#end if
			#end if
			if append_to_out:
				out.append((quadruplex_as_list_island, score[0]))
			#end if
		#end if
		score[0] -= results.pop().get_score()
	#end def __routine_get_paths_score

#end def g_graph


########################################
####		FUNCTIONS               ####
########################################

def is_island(seq, base, strt_idx, min_len, max_gaps_len=0, max_gaps=0):
	''' function to identify base island of lenght min_len,
	return bool value and information on the island '''
	if max_gaps_len and not max_gaps:
		max_gaps = min_len #maximum number of gaps possible
	#end if
	base_num, gaps_num, gaps_len, score = 1, 0, 0, 6
	last_base, increment, idx = base, 1, 1
	while (gaps_num <= max_gaps) and (gaps_len <= max_gaps_len):
		try:
			next_base = seq[strt_idx + idx]
			if base == next_base:
				base_num += 1
				score += 6
			else:
				gaps_len += 1
				if base == last_base:
					gaps_num += 1
					score -= 6
				#end if
				score -= 2
			#end if
			if base_num == min_len:
				return (1, base_num, gaps_num, gaps_len, score, strt_idx + idx)
			#end if
			last_base = next_base
			idx += 1
		except Exception:
			return (0, 0, 0, 0, 0, 0)
		#end try
	#end while
	return (0, 0, 0, 0, 0, 0)
#end def is_island

#### Island Scan ####
def island_scan(seq, base, min_len, max_gaps_len, max_gaps):
	''' function to scan a sequence seq searching
	for all the islands of length min_len, save into a list '''
	list_island = []
	name_i = 1
	for i, base_i in enumerate(seq):
		if base == base_i:
			check, base_num, gaps_num, gaps_len, score, end_idx = is_island(seq, base, i, min_len, max_gaps_len, max_gaps)
			if check:
				list_island.append(g_island(name_i, i, end_idx, base_num, gaps_num, gaps_len, score))
				name_i += 1
			#end if
		#end if
	#end for
	return list_island
#end def island_scan

def island_scan_range_run(seq, base, min_len, max_len, max_gaps_len, max_gaps, max_loop, of, length, all_flag, nocore_flag, min_perfect, faster_mode, normal_mode, symmetry, degenerancy, subst_matrix, indls, mismatch, bulge_only_mismatch):
	''' '''
	list_island, end_i, name_i, counter = [], 0, 1, 0
	for i, base_i in enumerate(seq):
		if base == base_i:
			list_island_i = routine_island_scan_range(seq, base, i, min_len, max_len, name_i, max_gaps_len, max_gaps, nocore_flag)
			# Run analysis if end block of possibly linked islands #
			if list_island_i:
				if (list_island_i[0].get_strt_idx() - end_i - 1) > max_loop:
					if counter >= length:
						graph_island = graph_build_island_length(list_island, max_loop, mismatch, bulge_only_mismatch)
						if verbose:
							sys.stderr.write('\tBuilt Graph up to {0}\n'.format(end_i))
							sys.stderr.flush()
						#end if
						if is_path(graph_island, list_island, length):
							if (faster_mode or length > 12) and not (normal_mode):
								search_length = 4
								# Calculate min_perfect #
								if min_perfect:
									min_perfect_per_four_islands = min_perfect / (length / 4)
									search_min_perfect = 4 if min_perfect_per_four_islands > 4 else min_perfect_per_four_islands
								else:
									search_min_perfect = min_perfect
								#end if
							else:
								search_length = length
								search_min_perfect = min_perfect
							#end if
							if not all_flag:
								print_paths_score(graph_island, list_island, seq, of, symmetry, degenerancy, subst_matrix, indls, mismatch, min_len, bulge_only_mismatch, search_length, search_min_perfect)
							else:
								print_paths(graph_island, list_island, seq, of, symmetry, degenerancy, subst_matrix, indls, mismatch, search_length)
							#end if
						#end if
					#end if
					# Reset variables #
					list_island, counter = [], 0
				#end if
				# Add new islands #
				counter += 1 #new islands group at index
				[list_island.append(island) for island in reversed(list_island_i)]
				end_i = list_island_i[-1].get_end_idx()
				name_i += 1
			#end if
		#end if
	#end for
	# Run analysis on last block of possibly linked islands #
	if list_island and counter >= length:
		graph_island = graph_build_island_length(list_island, max_loop, mismatch, bulge_only_mismatch)
		if verbose:
			sys.stderr.write('\tBuilt Graph up to {0}\n'.format(end_i))
			sys.stderr.flush()
		#end if
		if is_path(graph_island, list_island, length):
			if (faster_mode or length > 12) and not (normal_mode):
				search_length = 4
				# Calculate min_perfect #
				if min_perfect:
					min_perfect_per_four_islands = min_perfect / (length / 4)
					search_min_perfect = 4 if min_perfect_per_four_islands > 4 else min_perfect_per_four_islands
				else:
					search_min_perfect = min_perfect
				#end if
			else:
				search_length = length
				search_min_perfect = min_perfect
			#end if
			if not all_flag:
				print_paths_score(graph_island, list_island, seq, of, symmetry, degenerancy, subst_matrix, indls, mismatch, min_len, bulge_only_mismatch, search_length, search_min_perfect)
			else:
				print_paths(graph_island, list_island, seq, of, symmetry, degenerancy, subst_matrix, indls, mismatch, search_length)
			#end if
		#end if
	#end if
#end def island_scan_range_run

def routine_island_scan_range(seq, base, strt_idx, min_len, max_len, name_i, max_gaps_len=0, max_gaps=0, nocore_flag=False):
	''' function to identify base islands in range of length,
	return a list of island objects identified '''
	if max_gaps_len and not max_gaps:
		max_gaps = max_len #maximum number of gaps possible
	#end if
	out, base_num, gaps_num, gaps_len, score, island_len, core, perfect = [], 1, 0, 0, 6, min_len, False, True
	last_base, increment, idx = base, 1, 1
	while (gaps_num <= max_gaps) and (gaps_len <= max_gaps_len):
		try:
			next_base = seq[strt_idx + idx]
			if base == next_base:
				base_num += 1
				if base == last_base:
					core = True
				#end if
				score += 6
			else:
				gaps_len += 1
				perfect = False
				if base == last_base:
					gaps_num += 1
					score -= 6
				#end if
				score -= 2
			#end if
			if base_num == max_len:
				if not nocore_flag:
					if core:
						out.append(g_island(name_i, strt_idx, strt_idx + idx, base_num, gaps_num, gaps_len, score, perfect))
					#end if
				else:
					out.append(g_island(name_i, strt_idx, strt_idx + idx, base_num, gaps_num, gaps_len, score, perfect))
				#end if
				return out
			elif base_num == island_len:
				if not nocore_flag:
					if core:
						out.append(g_island(name_i, strt_idx, strt_idx + idx, base_num, gaps_num, gaps_len, score, perfect))
					#end if
				else:
					out.append(g_island(name_i, strt_idx, strt_idx + idx, base_num, gaps_num, gaps_len, score, perfect))
				#end if
				island_len +=  1
			#end if
			last_base = next_base
			idx += 1
		except Exception:
			return out
		#end try
	#end while
	return out
#end def routine_island_scan_range

#### Graph Build ####
def graph_build(list_island, max_loop):
	''' function to build a graph from the list of islands list_island,
	a max_loop is used to define which islands can be linked '''
	link_graph = graph(0)
	len_list_island = len(list_island)
	for i, island_i in enumerate(list_island):
		end_i, p = island_i.get_end_idx(), i + 1
		while p < len_list_island:
			island = list_island[p]
			start = island.get_strt_idx()
			if (start - end_i - 1) > max_loop:
				break
			elif (start > end_i):
				link_graph.add_node(island_i, island)
			#end if
			p += 1
		#end while
	#end for
	return link_graph
#end def graph_build

def graph_build_island_length(list_island, max_loop, mismatch, bulge_only_mismatch):
	''' function to build a graph from the list of islands list_island,
	a max_loop is used to define which islands can be linked '''
	link_graph = g_graph(0)
	len_list_island = len(list_island)
	for i, island_i in enumerate(list_island):
		end_i, p = island_i.get_end_idx(), i + 1
		while p < len_list_island:
			island = list_island[p]
			start = island.get_strt_idx()
			# Check mismatch #
			if mismatch and island_i.get_G_num() - 1 == island.get_G_num():
				if (start - end_i - 1) > max_loop:
					break
				elif (start > end_i):
					# Check if islands only has one internal mismatch #
					if bulge_only_mismatch:
						if island_i.get_gaps_len() == 0:
							link_graph.add_node(island_i, island)
						#end if
					else:
						if island.get_gaps_len() <= 1:
							link_graph.add_node(island_i, island)
						#end if
					#end if
				#end if
			elif mismatch and island_i.get_G_num() + 1 == island.get_G_num():
				if (start - end_i - 1) > max_loop:
					break
				elif (start > end_i):
					# Check if islands_i only has one internal mismatch #
					if bulge_only_mismatch:
						if island.get_gaps_len() == 0:
							link_graph.add_node(island_i, island)
						#end if
					else:
						if island_i.get_gaps_len() <= 1:
							link_graph.add_node(island_i, island)
						#end if
					#end if
				#end if
			else:
				if (start - end_i - 1) > max_loop:
					break
				elif (start > end_i) and (island_i.get_G_num() == island.get_G_num()):
					link_graph.add_node(island_i, island)
				#end if
			#end if
			p += 1
		#end while
	#end for
	return link_graph
#end def graph_build_island_length

#### Paths search ####
def print_paths(link_graph, list_island, seq, of, symmetry, degenerancy, subst_matrix, indls, mismatch, length=4):
	''' write all the possible quadruplex long length
	recorded in the link_graph to the output file of '''
	for island in list_island:
		for quadruplex_as_list_island in link_graph.get_paths(island, length):
			if symmetry:
				print_quadruplex_symmetry(quadruplex_as_list_island, seq, of, symmetry, degenerancy, subst_matrix, indls, mismatch)
			else:
				print_quadruplex(quadruplex_as_list_island, seq, of, mismatch)
			#end if
		#end for
	#end for
#end def print_paths

def is_path(link_graph, list_island, length=4):
	''' check if there is a path long length in the graph '''
	list_island_len = len(list_island)
	for i, island in enumerate(list_island):
		if verbose:
			sys.stderr.write('\r\t-> Check longest path [{:.0f}%]'.format(float(i)/list_island_len*100))
			sys.stderr.flush()
		#end if
		if link_graph.path_check(island, length):
			if verbose:
				sys.stderr.write('\r\t-> Check longest path [{0}%] - [OK]\n'.format(100))
				sys.stderr.flush()
			#end if
			return True
		#end if
	#end for
	if verbose:
		sys.stderr.write('\r\t-> Check longest path [{0}%] - [NOT FOUND]\n'.format(100))
		sys.stderr.flush()
	#end if
	return False
#end def is_path

def print_paths_score(link_graph, list_island, seq, of, symmetry, degenerancy, subst_matrix, indls, mismatch, min_len, bulge_only_mismatch, length=4, min_perfect=1):
	''' write all the best non contained quadruplex long length for each index
	recorded in the link_graph to the output file of '''
	last_name, last_score, printed, last_end_idx, end_idx_block = '', 0, False, 0, 0
	list_island_len = len(list_island)
	#print [(island.get_island(seq), island.get_strt_idx()) for island in list_island]
	for i, island in enumerate(list_island):
		if verbose:
			sys.stderr.write('\r\t-> Graph Navigation [{:.0f}%]'.format(float(i)/list_island_len*100))
			sys.stderr.flush()
		#end if
		current_name = island.get_name()
		# Check islands block per index #
		if last_name != current_name: #new block start
			last_name = current_name
			printed = False
		elif last_name == current_name and False == printed:
			pass
		else:
			continue
		#end if
		# Calculate and write quadruplex #
		lista_quadruplex = sorted(link_graph.get_paths_score(island, length, mismatch, min_len, bulge_only_mismatch, min_perfect), key=lambda(x,y): y, reverse=True)
		try:
			highest_score, end_idx_block = lista_quadruplex[0][1], 0
			for quadruplex_as_list_island, score in lista_quadruplex:
				if (highest_score == score):
					end_idx = quadruplex_as_list_island[-1].get_end_idx()
					if end_idx > end_idx_block:
						end_idx_block = end_idx
					#end if
					if (end_idx > last_end_idx) or (score >= last_score):
#					if (end_idx > last_end_idx) or (score > last_score):
						if symmetry:
							print_quadruplex_score_symmetry(quadruplex_as_list_island, seq, of, score, symmetry, degenerancy, subst_matrix, indls, mismatch)
						else:
							print_quadruplex_score(quadruplex_as_list_island, seq, of, score, mismatch)
						#end if
						printed = True
					#end if
				else:
					break
				#end if
			#end for
			if printed:
				last_end_idx = end_idx_block
				last_score = highest_score
			#end if
		except Exception:
			continue
		#end try
	#end for
	if verbose:
		sys.stderr.write('\r\t-> Graph Navigation [{0}%]\n'.format(100))
		sys.stderr.flush()
	#end if
#end def print_paths_score

#### Printing functions ####
def print_quadruplex(quadruplex_as_list_island, seq, of, mismatch):
	''' write a quadruplex on the output file of ''' #quadruplex, start_idx, end_idx, score
	#Print quadruplex
	score_total, idx = quadruplex_as_list_island[0].get_score(), quadruplex_as_list_island[0].get_end_idx() + 1
	of.write(quadruplex_as_list_island[0].get_island(seq)) #writing first island
	of.write('-')
	for island in quadruplex_as_list_island[1:-1]: #looping other island
		of.write(seq[idx:island.get_strt_idx()].lower())
		of.write('-')
		of.write(island.get_island(seq))
		of.write('-')
		score_total += island.get_score()
		idx = island.get_end_idx() + 1
	#end for
	of.write(seq[idx:quadruplex_as_list_island[-1].get_strt_idx()].lower())
	of.write('-')
	of.write(quadruplex_as_list_island[-1].get_island(seq)) #writing last island
	score_total += quadruplex_as_list_island[-1].get_score()
	#Print other information
	start_idx, end_idx = quadruplex_as_list_island[0].get_strt_idx(), quadruplex_as_list_island[-1].get_end_idx()
	if mismatch:
		of.write('\t{0}\t{1}\t{2}\t{3}\n'.format(score_total, start_idx, end_idx, max([island.get_G_num() for island in quadruplex_as_list_island])))
	else:
		of.write('\t{0}\t{1}\t{2}\t{3}\n'.format(score_total, start_idx, end_idx, island.get_G_num()))
	#end if
#end def print_quadruplex

def print_quadruplex_symmetry(quadruplex_as_list_island, seq, of, symmetry, degenerancy, subst_matrix, indls, mismatch):
	''' write a quadruplex on the output file of ''' #quadruplex, start_idx, end_idx, score
	#Print quadruplex
	list_alignments = []
	score_total, idx = quadruplex_as_list_island[0].get_score(), quadruplex_as_list_island[0].get_end_idx() + 1
	of.write(quadruplex_as_list_island[0].get_island(seq)) #writing first island
	of.write('-')
	for island in quadruplex_as_list_island[1:-1]: #looping other island
		sequence = seq[idx:island.get_strt_idx()].lower()
		of.write(sequence)
		#Check symmetry and adding score
		if len(sequence) > 6:
			max_score, alignment = check_symmetry(sequence, symmetry, degenerancy, subst_matrix, indls)
		else:
			max_score, alignment = 7 - len(sequence), ''
		#end if
		score_total += max_score
		if alignment:
			list_alignments.append(alignment)
		else:
			list_alignments.append('NA')
		#end if
		of.write('-')
		of.write(island.get_island(seq))
		of.write('-')
		score_total += island.get_score()
		idx = island.get_end_idx() + 1
	#end for
	sequence = seq[idx:quadruplex_as_list_island[-1].get_strt_idx()].lower()
	of.write(sequence)
	#Check symmetry and adding score
	if len(sequence) > 6:
		max_score, alignment = check_symmetry(sequence, symmetry, degenerancy, subst_matrix, indls)
	else:
		max_score, alignment = 7 - len(sequence), ''
	#end if
	score_total += max_score
	if alignment:
		list_alignments.append(alignment)
	else:
		list_alignments.append('NA')
	#end if
	of.write('-')
	of.write(quadruplex_as_list_island[-1].get_island(seq)) #writing last island
	score_total += quadruplex_as_list_island[-1].get_score()
	#Print other information
	start_idx, end_idx = quadruplex_as_list_island[0].get_strt_idx(), quadruplex_as_list_island[-1].get_end_idx()
	if mismatch:
		of.write('\t{0}\t{1}\t{2}\t{3}\t{4}\n'.format(score_total, start_idx, end_idx, max([island.get_G_num() for island in quadruplex_as_list_island]), ';'.join(list_alignments)))
	else:
		of.write('\t{0}\t{1}\t{2}\t{3}\t{4}\n'.format(score_total, start_idx, end_idx, island.get_G_num(), ';'.join(list_alignments)))
	#end if
#end def print_quadruplex_symmetry

def print_quadruplex_score(quadruplex_as_list_island, seq, of, score, mismatch):
	''' write a quadruplex on the output file of ''' #quadruplex, start_idx, end_idx, score
	#Print quadruplex
	idx = quadruplex_as_list_island[0].get_end_idx() + 1
	of.write(quadruplex_as_list_island[0].get_island(seq)) #writing first island
	of.write('-')
	for island in quadruplex_as_list_island[1:-1]: #looping other island
		of.write(seq[idx:island.get_strt_idx()].lower())
		of.write('-')
		of.write(island.get_island(seq))
		of.write('-')
		idx = island.get_end_idx() + 1
	#end for
	of.write(seq[idx:quadruplex_as_list_island[-1].get_strt_idx()].lower())
	of.write('-')
	of.write(quadruplex_as_list_island[-1].get_island(seq)) #writing last island
	#Print other information
	start_idx, end_idx = quadruplex_as_list_island[0].get_strt_idx(), quadruplex_as_list_island[-1].get_end_idx()
	if mismatch:
		of.write('\t{0}\t{1}\t{2}\t{3}\n'.format(score, start_idx, end_idx, max([island.get_G_num() for island in quadruplex_as_list_island])))
	else:
		of.write('\t{0}\t{1}\t{2}\t{3}\n'.format(score, start_idx, end_idx, island.get_G_num()))
	#end if
#end def print_quadruplex_score

def print_quadruplex_score_symmetry(quadruplex_as_list_island, seq, of, score, symmetry, degenerancy, subst_matrix, indls, mismatch):
	''' write a quadruplex on the output file of ''' #quadruplex, start_idx, end_idx, score
	#Print quadruplex
	list_alignments = []
	idx = quadruplex_as_list_island[0].get_end_idx() + 1
	of.write(quadruplex_as_list_island[0].get_island(seq)) #writing first island
	of.write('-')
	for island in quadruplex_as_list_island[1:-1]: #looping other island
		sequence = seq[idx:island.get_strt_idx()].lower()
		of.write(sequence)
		#Check symmetry and adding score
		if len(sequence) > 6:
			max_score, alignment = check_symmetry(sequence, symmetry, degenerancy, subst_matrix, indls)
		else:
			max_score, alignment = 7 - len(sequence), ''
		#end if
		score += max_score
		if alignment:
			list_alignments.append(alignment)
		else:
			list_alignments.append('NA')
		#end if
		of.write('-')
		of.write(island.get_island(seq))
		of.write('-')
		idx = island.get_end_idx() + 1
	#end for
	sequence = seq[idx:quadruplex_as_list_island[-1].get_strt_idx()].lower()
	of.write(sequence)
	#Check symmetry and adding score
	if len(sequence) > 6:
		max_score, alignment = check_symmetry(sequence, symmetry, degenerancy, subst_matrix, indls)
	else:
		max_score, alignment = 7 - len(sequence), ''
	#end if
	score += max_score
	if alignment:
		list_alignments.append(alignment)
	else:
		list_alignments.append('NA')
	#end if
	of.write('-')
	of.write(quadruplex_as_list_island[-1].get_island(seq)) #writing last island
	#Print other information
	start_idx, end_idx = quadruplex_as_list_island[0].get_strt_idx(), quadruplex_as_list_island[-1].get_end_idx()
	if mismatch:
		of.write('\t{0}\t{1}\t{2}\t{3}\t{4}\n'.format(score, start_idx, end_idx, max([island.get_G_num() for island in quadruplex_as_list_island]), ';'.join(list_alignments)))
	else:
		of.write('\t{0}\t{1}\t{2}\t{3}\t{4}\n'.format(score, start_idx, end_idx, island.get_G_num(), ';'.join(list_alignments)))
	#end if
#end def print_quadruplex_score_symmetry

#### Alignment ####
def check_symmetry(sequence, symmetry, degenerancy, matrix, indls):
	''' modified Needleman-Wunsch algorithm '''
	k = len(sequence)
#	indls = -1
#	max_errors = k * degenerancy // 100
	max_i, max_j, max_score = 0, 0, 0
	seq, alignment = sequence.upper().replace('U', 'T'), ''
	#Check Iupac not supported
	if ('N' in seq or 'R' in seq or 'Y' in seq or 'S' in seq or 'W' in seq or 'K' in seq \
		or 'D' in seq or 'M' in seq or 'H' in seq or 'B' in seq or 'V' in seq or 'I' in seq):
		return 0, ''
	#end if
	score_matrix = np.zeros(shape=(k + 1, k + 1), dtype=np.int)
	subst_idx = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
	complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
	if symmetry == 'mirror':
		subst_matrix = np.matrix([[1, -1, -1, -1], [-1, 1, -1, -1], [-1, -1, 1, -1], [-1, -1, -1, 1]])
	elif symmetry == 'palindrome':
		subst_matrix = np.matrix([[-1, -1, -1, 1], [-1, -1, 1, -1], [-1, 1, -1, -1], [1, -1, -1, -1]])
	elif symmetry == 'mixed':
		#subst_matrix = np.matrix([[1, -1, -1, 1], [-1, 1, 1, -1], [-1, 1, 1, -1], [1, -1, -1, 1]])
		subst_matrix = np.matrix([[1, -2, -2, 2], [-2, 1, 2, -2], [-2, 2, 1, -2], [2, -2, -2, 1]])
	else:
		subst_matrix = matrix
	#end if
	#Initialize score_matrix
	for i in range(k + 1):
		score_matrix[0, i] = i * indls
		score_matrix[i, 0] = i * indls
	#end for
	#Compute matrix
	max_score_def = False
	for i in range(1, k + 1):
		for j in range(1, k + 1):
			if (i + j) > k:
				break
			#end if
			score_up = score_matrix[i - 1, j] + indls
			score_left = score_matrix[i, j - 1] + indls
			score_diag = score_matrix[i - 1, j - 1] + subst_matrix[subst_idx[seq[i - 1]], subst_idx[seq[k - j]]]
			score_matrix[i, j] = max(score_up, score_left, score_diag)
			#Track highest score
			if k == (i + j): #along the diagonal
				if not max_score_def:
					max_score, max_i, max_j = score_matrix[i, j], i, j
					max_score_def = True
				elif max_score < score_matrix[i, j]:
					max_score, max_i, max_j = score_matrix[i, j], i, j
				#end if
			#end if
		#end for
	#end for
	#Traceback
	i, j = max_i, max_j
	while (i > 0) or (j > 0):
		if (i > 0) and (j > 0) and (score_matrix[i, j] == score_matrix[i - 1, j - 1] + subst_matrix[subst_idx[seq[i - 1]], subst_idx[seq[k - j]]]):
			if seq[i - 1] == complement[seq[k - j]]: #Watson-Crick
				alignment += 'W'
			elif seq[i - 1] == seq[k - j]: #Hoogsteen
				alignment += 'H'
			else:
				alignment += 'm'
			#end if
			i -= 1
			j -= 1
		elif (i > 0) and (score_matrix[i, j] == score_matrix[i - 1, j] + indls):
			alignment += 'u'
			i -= 1
		else:
			alignment += 'l'
			j -= 1
		#end if
	#end while
	penalty, edited_alignment = edit_alignment(seq, alignment, subst_matrix, subst_idx, indls)
	return max_score + penalty, edited_alignment
#end def check_symmetry

def edit_alignment(seq, alignment, subst_matrix, subst_idx, indls):
	''' '''
	penalty, edited_alignment, k, l = 0, '', len(alignment), alignment.count('l')
	j = l
	if alignment[0] == 'l': l -= 1
	if alignment[0] == 'u': j += 1
	for i in range(2):
		if alignment[i] == 'l' or alignment[i] == 'u':
			penalty -= indls
			edited_alignment += alignment[i]
		elif alignment[i] == 'W':
			penalty -= subst_matrix[subst_idx[seq[k - l - i - 1]], subst_idx[seq[k - j + i]]]
			edited_alignment += 'm'
		elif alignment[i] == 'H':
			penalty -= subst_matrix[subst_idx[seq[k - l - i - 1]], subst_idx[seq[k - j + i]]]
			edited_alignment += 'm'
		else:
			penalty -= subst_matrix[subst_idx[seq[k - l - i - 1]], subst_idx[seq[k - j + i]]]
			edited_alignment += alignment[i]
		#end if
	#end for
	if k > 2: edited_alignment += alignment[2:]
	return penalty, edited_alignment
#end def edit_alignment


########################################
####		MAIN                    ####
########################################

def main(args): # use as args['name']

	## Variables ##
	min_len = int(args['minlen']) if args['minlen'] else 2
	max_len = int(args['maxlen']) if (args['maxlen']) and (args['maxlen'] > min_len) else min_len
	max_loop = int(args['maxloop']) if args['maxloop'] else 12
	max_gaps = int(args['gapnum']) if args['gapnum'] else 0
	max_gaps_len = int(args['gaplen']) if args['gaplen'] else 0
	base = args['base'].upper() if args['base'] else 'G'
	island_num = int(args['islandnum']) if args['islandnum'] else 4
	base_counter, checked_min_bases = 0, False
	nocore = args['nocore']
	allresult = False
	# Check gap #
	if max_gaps and not max_gaps_len:
		sys.exit('\ninput error: only gaps number (-g) have been selected, please select a length for the gaps (-l)\n')
	#end if
	# Check Perfect #
	if args['noperfect']:
		min_perfect = 0
	else:
		if args['perfect'] and int(args['perfect']) > island_num:
			sys.exit('\ninput error: perfect islands required are higher than the number of islands required\n')
		#end if
		min_perfect = int(args['perfect']) if args['perfect'] else 1
	#end if
	# Check faster mode #
	if args['fastermode']:
		if island_num <= 4:
			sys.exit('\ninput error: faster mode can be only applied to the search of more than four islands\n')
		#end if
	#end if
	# Check symmetry mode #
	subst_matrix = np.zeros(shape=(4, 4), dtype=np.int)
	indls = -1
	if (args['simmetrymirror'] and args['simmetrypalindrome']) or (args['simmetrymirror'] and args['simmetrymixed']) or (args['simmetrypalindrome'] and args['simmetrymixed']):
		sys.exit('\ninput error: multiple symmetry modes selected, please select only one at the time\n')
	elif args['simmetrymirror'] or args['simmetrypalindrome'] or args['simmetrymixed']:
		if args['simmetrymirror']:
			symmetry, degenerancy = 'mirror', 100
		elif args['simmetrypalindrome']:
			symmetry, degenerancy = 'palindrome', 100
		else:
			symmetry, degenerancy, indls = 'mixed', 100, -3
		#end if
	else:
		symmetry, degenerancy = False, 0
	#end if
	# Check custom matrix #
	if args['simmetrycustom']:
		if symmetry:
			sys.exit('\ninput error: multiple symmetry modes selected, please select only one at the time\n')
		else:
			subst_idx = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
			symmetry, degenerancy = 'custom', 100
			read = open(args['simmetrycustom'])
			for line in read:
				par, score = line.rstrip().split('=')
				if par.strip() == 'gap_open':
					indls = int(score.strip())
				else:
					subst_matrix[subst_idx[par.strip()[0]], subst_idx[par.strip()[1]]] = int(score.strip())
					subst_matrix[subst_idx[par.strip()[1]], subst_idx[par.strip()[0]]] = int(score.strip())
				#end if
			#end for
			read.close()
		#end if
	#end if
	# Override parameters when mismatch selected #
	bulge_only_mismatch = False
	if args['mismatch']:
		if min_len < 3:
			sys.exit('\ninput error: mismatched-islands can be only evaluated for a minimum islands length of at least three\n')
		elif int(args['mismatch']) > (island_num - 1):
			sys.exit('\ninput error: number of mismatched-islands allowed is higher than the number of islands required\n')
		else:
			min_len -= 1
			nocore = True
			mismatch = int(args['mismatch'])
			if not max_gaps_len:
				max_gaps_len = 1
				bulge_only_mismatch = True
			#end if
		#end if
	else:
		mismatch = 0
	#end if

	## Opening output file ##
	fo = open(args['output'], 'w')
	if symmetry:
		fo.write('#quadruplex\tscore\tstart\tend\tisland_len\tloops_symmetry\n')
	else:
		fo.write('#quadruplex\tscore\tstart\tend\tisland_len\n')
	#end if

	## Reading the input sequences and printing quadruplex ##
	with open(args['inputsequence'], 'r') as fi:
		first_header = True
		for line in fi:
			if line.startswith('>') and first_header:
				first_header = False
				header, seq = line.rstrip()[1:], ''
			elif line.startswith('>') and not first_header:
				# Printing header #
				sys.stderr.write('Analyzing sequence - {0}\n'.format(header))
				sys.stderr.flush()
				fo.write('>{0}\n'.format(header))

				# Analyzing the sequence #
				if False:
				#if not checked_min_bases:
					island_scan_range_run(seq, base, 2, max_len, max_gaps_len, max_gaps, max_loop, fo, island_num,
										allresult, nocore, min_perfect, args['fastermode'], args['normalmode'],
										symmetry, degenerancy, subst_matrix, indls, mismatch, bulge_only_mismatch)
				else:
					island_scan_range_run(seq, base, min_len, max_len, max_gaps_len, max_gaps, max_loop, fo, island_num,
										allresult, nocore, min_perfect, args['fastermode'], args['normalmode'],
										symmetry, degenerancy, subst_matrix, indls, mismatch, bulge_only_mismatch)
				#end if

				# Reset Variables #
				header, seq = line.rstrip()[1:], ''
				base_counter, checked_min_bases = 0, False
			elif checked_min_bases:
				seq += line.rstrip().upper()
			else:
				for line_i in line.rstrip().upper():
					if base == line_i:
						base_counter += 1
					#end if
					if base_counter >= (min_len * island_num):
						checked_min_bases = True
					#end if
				#end for
				seq += line.rstrip().upper()
			#end if
		#end for
	#end with

	## Analyzing and printing last sequence ##
	sys.stderr.write('Analyzing LAST sequence - {0}\n'.format(header))
	sys.stderr.flush()
	fo.write('>{0}\n'.format(header))

	if False:
	#if not checked_min_bases:
		island_scan_range_run(seq, base, 2, max_len, max_gaps_len, max_gaps, max_loop, fo, island_num,
							allresult, nocore, min_perfect, args['fastermode'], args['normalmode'],
							symmetry, degenerancy, subst_matrix, indls, mismatch, bulge_only_mismatch)
	else:
		island_scan_range_run(seq, base, min_len, max_len, max_gaps_len, max_gaps, max_loop, fo, island_num,
							allresult, nocore, min_perfect, args['fastermode'], args['normalmode'],
							symmetry, degenerancy, subst_matrix, indls, mismatch, bulge_only_mismatch)
	#end if

	## Closing output file ##
	fo.close()

# end def main


if __name__ == '__main__':

	parser = argparse.ArgumentParser(description='Welcome to QPARSE.\nCopyright(C) 2019 Michele Berselli\nberselli.michele@gmail.com',
									formatter_class=argparse.RawTextHelpFormatter)

	parser.add_argument('-i','--inputsequence', help='input sequence as fasta/multifasta', required=True)
	parser.add_argument('-o','--output', help='output file', required=True)
	parser.add_argument('-m','--minlen', help='minimum islands length [2]', required=False)
	parser.add_argument('-M','--maxlen', help='maximum islands length [2]', required=False)
	parser.add_argument('-L','--maxloop', help='maximum loops length [12]', required=False)
	parser.add_argument('-g','--gapnum', help='maximum number of gaps per island [0]', required=False)
	parser.add_argument('-l','--gaplen', help='maximum length of gaps per island [0]', required=False)
	parser.add_argument('-b','--base', help='base to use for search [G]', required=False)
	parser.add_argument('-n','--islandnum', help='number of consecutive islands [4]', required=False)
	parser.add_argument('-p','--perfect', help='minimum number of perfect islands (no-gaps) required [1]', required=False)
	parser.add_argument('-x','--mismatch', help='evaluate also islands containing a mismatch,\nset the maximum number of mismatched-islands allowed [0]', required=False)
#	parser.add_argument('-all','--allresult', help='show all possible patterns, also overlapping and suboptimal\n[caution: the output can be huge]', action='store_true', required=False)
	parser.add_argument('-nocore','--nocore', help='detect also islands without at least two consecutive bases (e.g. GaGcG)', action='store_true', required=False)
	parser.add_argument('-noperfect','--noperfect', help='detect also patterns without any perfect island', action='store_true', required=False)
	parser.add_argument('-fast','--fastermode',
						help='faster search mode that can be applied for islandnum > 4\n[only regions with at least --islandnum islands are evaluated\nto build the graph, however, the graph is navigated\nsearching for patterns of four islands to reduce the\ncombinatorial complexity and speed up the analysis]',
						action='store_true', required=False)
	parser.add_argument('-normal','--normalmode',
						help='when searching for more than 12 islands --fastermode is used by default,\nthis flag override --fastermode and the graph is navigated\nsearching for patterns of --islandnum islands\n[caution: the analysis can be computationally expensive]',
						action='store_true', required=False)
	parser.add_argument('-sM','--simmetrymirror',
						help='evaluate the symmetry of the long loops (>= 7 bp) to improve the score,\nMIRROR symmetry is considered\n[allows to detect longer loops with mirror properties that can form Hoogsteen-hairpins]',
						action='store_true', required=False)
	parser.add_argument('-sP','--simmetrypalindrome',
						help='evaluate the symmetry of the long loops (>= 7 bp) to improve the score,\nPALINDROMIC symmetry is considered\n[allows to detect longer loops with palindromic properties that can form canonical-hairpins]',
						action='store_true', required=False)
	parser.add_argument('-sX','--simmetrymixed',
						help='evaluate the symmetry of the long loops (>= 7 bp) to improve the score,\nMIXED MIRROR-PALINDROMIC symmetry is considered\n[allows to detect longer loops with mirror and palindromic properties that can form mixed-hairpins]',
						action='store_true', required=False)
	parser.add_argument('-sC','--simmetrycustom',
						help='evaluate the symmetry of the long loops (>= 7 bp) to improve the score,\nthe custom substitution-matrix specified is used\n[allows to detect longer loops with symmetrical properties that can form hairpins]',
						required=False)
	parser.add_argument('-v','--verbose',
						help='activate verbose mode to show graph building status',
						action='store_true', required=False)

	args = vars(parser.parse_args())

	verbose = args['verbose']

	main(args)

# end if
