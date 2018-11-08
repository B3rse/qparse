#!/usr/bin/env python


######################################################################################
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


# Import libraries
import sys, argparse


# Functions definition
def routine_print_gff(of, seqID, seq, start, end, island_len, score='.'):
	''' '''
	of.write('{0}\tQPARSE\t{1}\t{2}\t{3}\t{4}\t{5}\t.\tID={6};island={7}\n'.format(
																			seqID,
																			'quadruplex',
																			start + 1,
																			end + 1,
																			score,
																			'+',
																			seq,
																			island_len
																				))
#end def routine_print_gff

def routine_print_tsv(of, seqID, seq, start, end, island_len, score='.'):
	''' '''
	of.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n'.format(seqID, start, end, island_len, score, seq)) #seqID	start	end	island_len	score	quadruplex
#end def routine_print_tsv

def merge(fi, fo, gff, maxLoop):
	''' '''
	if not gff:
		fo.write('#seqID\tstart\tend\tisland_len\tscore\tquadruplex\n')
	#end if
	first_ID = True
	for line in fi:
		if not line.startswith('#'):
			if line.startswith('>'):
				if first_ID:
					first_ID = False
				else:
					if gff and not first:
						routine_print_gff(fo, seqID, seq, start_i, end_i, island_len_i, score_i)
					elif not first:
						routine_print_tsv(fo, seqID, seq, start_i, end_i, island_len_i, score_i)
					#end if
				#end if
				seqID = line.rstrip()[1:]
				first = True
			else:
				line_splitted = line.rstrip().split()
				g4ID, start, end, island_len, score = line_splitted[0], int(line_splitted[2]), int(line_splitted[3]), int(line_splitted[4]), int(line_splitted[1])
				check = True
				if maxLoop:
					check = check_maxLoop(g4ID, maxLoop)
				#end if
				if first and check:
					first = False
					start_i, end_i, seq, island_len_i, score_i = start, end, '', island_len, score
					seq += g4ID.replace('-', '').upper()
				elif check:
					if start > end_i:
						if gff:
							routine_print_gff(fo, seqID, seq, start_i, end_i, island_len_i, score_i)
						else:
							routine_print_tsv(fo, seqID, seq, start_i, end_i, island_len_i, score_i)
						#end if
						start_i, end_i, seq, island_len_i, score_i = start, end, '', island_len, score
						seq += g4ID.replace('-', '').upper()
					elif end > end_i:
						idx = end - end_i
						end_i = end
						seq += g4ID.replace('-', '').upper()[-idx:]
					#end if
					if island_len > island_len_i:
						island_len_i = island_len
					#end if
					if score > score_i:
						score_i = score
					#end if
				#end if
			#end if
		#end if
	#end for
	if gff and not first:
		routine_print_gff(fo, seqID, seq, start_i, end_i, island_len_i, score_i)
	elif not first:
		routine_print_tsv(fo, seqID, seq, start_i, end_i, island_len_i, score_i)
	#end if
#end def merge

def max_number(fi, fo, gff, maxLoop):
	''' '''
	if not gff:
		fo.write('#seqID\tstart\tend\tisland_len\tscore\tquadruplex\n')
	#end if
	first_ID = True
	for line in fi:
		if not line.startswith('#'):
			if line.startswith('>'):
				if first_ID:
					first_ID = False
				else:
					if gff and not printed and not first:
						routine_print_gff(fo, seqID, seq_i, start_i, end_i, island_len_i, score_i)
						printed = True
					elif not printed and not first:
						routine_print_tsv(fo, seqID, seq_i, start_i, end_i, island_len_i, score_i)
						printed = True
					#end if
				#end if
				seqID = line.rstrip()[1:]
				first, printed = True, True
			else:
				line_splitted = line.rstrip().split()
				g4ID, start, end, island_len, score = line_splitted[0], int(line_splitted[2]), int(line_splitted[3]), int(line_splitted[4]), int(line_splitted[1])
				check = True
				if maxLoop:
					check = check_maxLoop(g4ID, maxLoop)
				#end if
				if first and check:
					first = False
					start_i, end_i, seq_i, island_len_i, score_i, printed = start, end, g4ID, island_len, score, False
				elif check:
					if start == start_i:
						if end < end_i:
							start_i, end_i, seq_i, island_len_i, score_i = start, end, g4ID, island_len, score
						#end if
					else:
						if gff and not printed:
							routine_print_gff(fo, seqID, seq_i, start_i, end_i, island_len_i, score_i)
							printed = True
						elif not printed:
							routine_print_tsv(fo, seqID, seq_i, start_i, end_i, island_len_i, score_i)
							printed = True
						#end if
						if start > end_i:
							start_i, end_i, seq_i, island_len_i, score_i, printed = start, end, g4ID, island_len, score, False
						#end if
					#end if
				#end if
			#end if
		#end if
	#end for
	if gff and not printed and not first:
		routine_print_gff(fo, seqID, seq_i, start_i, end_i, island_len_i, score_i)
		printed = True
	elif not printed and not first:
		routine_print_tsv(fo, seqID, seq_i, start_i, end_i, island_len_i, score_i)
		printed = True
	#end if
#end def max

def normal(fi, fo, gff, maxLoop):
	''' '''
	if not gff:
		fo.write('#seqID\tstart\tend\tisland_len\tscore\tquadruplex\n')
	#end if
	for line in fi:
		if not line.startswith('#'):
			if line.startswith('>'):
				seqID = line.rstrip().split()[0][1:]
			else:
				line_splitted = line.rstrip().split()
				seq, score, start, end, island_len = line_splitted[0], int(line_splitted[1]), int(line_splitted[2]), int(line_splitted[3]), int(line_splitted[4])
				printa = True
				if maxLoop:
					printa = check_maxLoop(seq, maxLoop)
				#end if
				if gff and printa:
					routine_print_gff(fo, seqID, seq, start, end, island_len, score)
				elif printa:
					routine_print_tsv(fo, seqID, seq, start, end, island_len, score)
				#end if
			#end if
		#end if
	#end for
#end def normal

def score(fi, fo, alignment, maxLoop):
	''' '''
	fo.write('#seqID\tstart\tend\tisland_len\tscore\tquadruplex\n')
	dict_quadruplex, dict_quadruplex_aln = {}, {}
	c = 0
	for line in fi:
		if not line.startswith('#'):
			if line.startswith('>'):
				c += 1
				seqID = line.rstrip().split()[0][1:]
			else:
				line_splitted = line.rstrip().split()
				printa = True
				if len(line_splitted) == 6 and alignment:
					seq, score, start, end, island_len, aln = line_splitted[0], int(line_splitted[1]), int(line_splitted[2]), int(line_splitted[3]), int(line_splitted[4]), line_splitted[5]
					if maxLoop:
						printa = check_maxLoop(seq, maxLoop)
					#end if
					if printa:
						dict_quadruplex_aln.setdefault((seqID, c), {})
						dict_quadruplex_aln[(seqID, c)].setdefault((seq, start), [score, end, island_len, aln])
					#end if
				elif len(line_splitted) < 6 and alignment:
					sys.exit('\nruntime error: input file was generated without checking for loop symmetry\n')
				else:
					seq, score, start, end, island_len = line_splitted[0], int(line_splitted[1]), int(line_splitted[2]), int(line_splitted[3]), int(line_splitted[4])
					if maxLoop:
						printa = check_maxLoop(seq, maxLoop)
					#end if
					if printa:
						dict_quadruplex.setdefault((seqID, c), {})
						dict_quadruplex[(seqID, c)].setdefault((seq, start), [score, end, island_len])
					#end if
				#end if
			#end if
		#end if
	#end for

	if dict_quadruplex:
		for (seqID, c), dict_seqID in sorted(dict_quadruplex.iteritems(), key=lambda (x, y): x[1]):
			for (seq, start), (score, end, island_len) in sorted(dict_seqID.iteritems(), key=lambda (x, y): y[0], reverse = True):
				routine_print_tsv(fo, seqID, seq, start, end, island_len, score)
			#end for
		#end for
	else:
		for (seqID, c), dict_seqID in sorted(dict_quadruplex_aln.iteritems(), key=lambda (x, y): x[1]):
			for (seq, start), (score, end, island_len, aln) in sorted(dict_seqID.iteritems(), key=lambda (x, y): y[0], reverse = True):
				routine_print_tsv(fo, seqID, seq, start, end, island_len, score)
				lista_aln = aln.split(';')
				lista_loop = seq.split('-') #[_, loop1, _, loop2, _, ...]
				i, i_loop = 1, 1
				for aln in lista_aln:
					fo.write('\tLoop_{0}:\n'.format(i))
					routine_print_align(fo, lista_loop[i_loop], aln[::-1])
					i += 1
					i_loop += 2
				#end for
			#end for
		#end for
	#end if
#end def score

def routine_print_align(of, seq, align):

	seq_len = len(seq)
	# Write first sequence
	of.write('\t\t')
	i = 0
	for l in align:
		if l == 'W' or l == 'm' or l == 'H' or l == 'u':
			of.write('{0}'.format(seq[i]))
			i += 1
		else:
			of.write('-')
		#end if
	#end for
	of.write('\n')

	# Write align
	of.write('\t\t')
	for l in align:
		if l == 'W':
			of.write('|')
		elif l == 'H':
			of.write(':')
		else:
			of.write(' ')
		#end if
	#end for
	of.write('\n')

	# Write second sequence
	of.write('\t\t')
	i = 0
	for l in align:
		if l == 'W' or l == 'm' or l == 'H' or l == 'l':
			of.write('{0}'.format(seq[seq_len - i - 1]))
			i += 1
		else:
			of.write('-')
		#end if
	#end for
	of.write('\n')
#end def routine_print_align

def check_maxLoop(seq, maxLoop):
	''' '''
	lista_loop = seq.split('-') #[_, loop1, _, loop2, _, ...]
	Loop_counter = 0
	for i, s_i in enumerate(lista_loop):
		if (i % 2) == 1:
			if len(s_i) > 5:
				Loop_counter += 1
			#end if
			if Loop_counter > maxLoop:
				return False
			#end if
		#end if
	#end for
	return True
#end def check_maxLoop

# Main
def main(args):

	## Opening files
	fi = open(args['inputfile'], 'r')
	fo = open(args['outputfile'], 'w')

	if args['alignment']:
		if args['gff'] or args['merge'] or args['max']:
			sys.exit('\ninput error: alignment mode is incompatible with other parsing modes, please select only one\n')
		#end if
	#end if
	
	if args['maxLoop']:
		maxLoop = int(args['maxLoop'])
	else:
		maxLoop = 0
	#end if

	## Parser
	if args['alignment'] or args['score']:
		score(fi, fo, args['alignment'], maxLoop)
	elif args['merge']:
		merge(fi, fo, args['gff'], maxLoop)
	elif args['max']:
		max_number(fi, fo, args['gff'], maxLoop)
	else:
		normal(fi, fo, args['gff'], maxLoop)
	#end if

	## Closing output file
	fo.close()
	fi.close()
#end def main


if __name__ == '__main__':

	parser = argparse.ArgumentParser(description='Script to format the output from QPARSE', 
									formatter_class=argparse.RawTextHelpFormatter)

	parser.add_argument('-i','--inputfile', help='output file from QPARSE as input', required=True)
	parser.add_argument('-o','--outputfile', help='file to store formatted output', required=True)
	parser.add_argument('-g','--gff', help='output to gff format instead of tsv', action='store_true', required=False)
	parser.add_argument('-m','--merge', help='merge overlapping patterns\n(return longest regions of consecutive islands whithin loop distance)', action='store_true', required=False)
	parser.add_argument('-x','--max', help='return the maximum number of non overlapping patterns', action='store_true', required=False)
	parser.add_argument('-s','--score', help='order results by score', action='store_true', required=False)
	parser.add_argument('-a','--alignment', help='order results by score,\nreturn the optimal alignment calculated for the linking loops', action='store_true', required=False)
	parser.add_argument('-maxLoop','--maxLoop', help='maximum number of long loops (>= 6 bp)', required=False)

	args = vars(parser.parse_args())

	main(args)
#end if

