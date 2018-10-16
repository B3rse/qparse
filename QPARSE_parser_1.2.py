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

def merge(fi, fo, gff):
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
						routine_print_gff(fo, seqID, seq, start_i, end_i, island_len_i)
					elif not first:
						routine_print_tsv(fo, seqID, seq, start_i, end_i, island_len_i)
					#end if
				#end if
				seqID = line.rstrip()[1:]
				first = True
			else:
				line_splitted = line.rstrip().split()
				g4ID, start, end, island_len = line_splitted[0], int(line_splitted[2]), int(line_splitted[3]), int(line_splitted[4])
				if first:
					first = False
					start_i, end_i, seq, island_len_i = start, end, '', island_len
					seq += g4ID.replace('-', '').upper()
				else:
					if start > end_i:
						if gff:
							routine_print_gff(fo, seqID, seq, start_i, end_i, island_len_i)
						else:
							routine_print_tsv(fo, seqID, seq, start_i, end_i, island_len_i)
						#end if
						start_i, end_i, seq, island_len_i = start, end, '', island_len
						seq += g4ID.replace('-', '').upper()
					elif end > end_i:
						idx = end - end_i
						end_i = end
						seq += g4ID.replace('-', '').upper()[-idx:]
						if island_len < island_len_i:
							island_len_i = island_len
						#end if
					#end if
				#end if
			#end if
		#end if
	#end for
	if gff and not first:
		routine_print_gff(fo, seqID, seq, start_i, end_i, island_len_i)
	elif not first:
		routine_print_tsv(fo, seqID, seq, start_i, end_i, island_len_i)
	#end if
#end def merge

def max_number(fi, fo, gff):
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
				if first:
					first = False
					start_i, end_i, seq_i, island_len_i, score_i, printed = start, end, g4ID, island_len, score, False
				else:
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

def normal(fi, fo, gff):
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
				if gff:
					routine_print_gff(fo, seqID, seq, start, end, island_len, score)
				else:
					routine_print_tsv(fo, seqID, seq, start, end, island_len, score)
				#end if
			#end if
		#end if
	#end for
#end def normal


# Main
def main(args):

	## Opening files
	fi = open(args['inputfile'], 'r')
	fo = open(args['outputfile'], 'w')

	## Parser
	if args['merge']:
		merge(fi, fo, args['gff'])
	elif args['max']:
		max_number(fi, fo, args['gff'])
	else:
		normal(fi, fo, args['gff'])
	#end if

	## Closing output file
	fo.close()
	fi.close()
#end def main


if __name__ == '__main__':

	parser = argparse.ArgumentParser(description='Script to format the output from QPARSE')

	parser.add_argument('-i','--inputfile', help='output file from QPARSE as input', required=True)
	parser.add_argument('-o','--outputfile', help='file to store formatted output', required=True)
	parser.add_argument('-g','--gff', help='output to gff format inseatd of tsv', action='store_true', required=False)
	parser.add_argument('-m','--merge', help='merge overlapping solutions', action='store_true', required=False)
	parser.add_argument('-x','--max', help='return the maximum number of non overlapping solutions', action='store_true', required=False)
	parser.add_argument('-a','--alignment', help='return the alignment calculated for the linking loops', action='store_true', required=False)

	args = vars(parser.parse_args())

	main(args)
#end if

