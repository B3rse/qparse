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
def routine_print_gff(seq_ID, g4_ID, score, start, end, island_len, of):

	of.write('{0}\tQUASAR\t{1}\t{2}\t{3}\t{4}\t{5}\t.\tID={6};island={7}\n'.format(
																			seq_ID, 
																			'quadruplex',
																			start + 1,
																			end + 1,
																			score,
																			'+',
																			g4_ID,
																			island_len
																				))
#end def routine_print_gff

def main(args):

	## Opening output file
	of = open(args['outputfile'], 'w')

	## Reading input file
	with open(args['inputfile']) as fi:
		for line in fi:
			if not line.startswith('#'):
				if line.startswith('>'):
					seq_ID = line.rstrip().split()[0][1:]
				else:
					line_splitted = line.rstrip().split()
					g4_ID, score, start, end, island_len = line_splitted[0], int(line_splitted[1]), int(line_splitted[2]), int(line_splitted[3]), int(line_splitted[4])
					routine_print_gff(seq_ID, g4_ID, score, start, end, island_len, of)
				#end if
			#end if
		#end for
	#end with

	## Closing output file
	of.close()
#end def main


if __name__ == '__main__':

	parser = argparse.ArgumentParser(description='Script to format the output from NeSSIe')

	parser.add_argument('-i','--inputfile', help='output file from NeSSIe as input', required=True)
	parser.add_argument('-o','--outputfile', help='file to store formatted output', required=True)

	args = vars(parser.parse_args())

	main(args)
#end if

