# QPARSE

**QPARSE** is a Python program (Python v. 2.7) that allows to search for *exact* or *degenerate* putative patterns associated with the formation of four-stranded non-canonical nucleic-acids structures (quadruplex or PQS) and multimeric quadruplex structures [[1](#references), [2](#references)].

The tool can assess the symmetrical properties of the linking loops to evaluate the potential of longer loops (>= 7 nt) to form stem-loop structures that stabilize the quadruplex [[2](#references), [3](#references)].

QPARSE exploits an exhaustive graphs-based algorithm that allows to model all the possible patterns considering all the combinations of islands that are detected based on the specified parameters.

### Features:

  - The tool can detect both *exact* or *degenerate* islands (i.e. islands containing any desired combination of bases and bulges, and mismatched islands).
  - The tool can detect longer PQS with more than four consecutive islands that can potentially form multimeric quadruplexes.
  - The tool is exhaustive in the search and all the possible combinations of the detected islands are considered.
  - The tool can assess the symmetrical properties (mirror, palindromic or a combination of both) of the linking loops to evaluate the potential formation of hairpins that stabilize longer loops (>= 7 nt). It is also possible to use a custom matrix for the calculation of the optimal self-alignment.
  - The search is not limited to G or C but can be extended to any other base.
  - A parser is provided to refine the raw output and convert it into a tsv or gff format. The parser also allows the visualization of the symmetrical properties of the loops in a blast-like format.

## **Contacts**
Michele Berselli, <berselli.michele@gmail.com>

Feel free to contact me for any inquiry, unexpected behaviour or support for your analyses.

A web-server is available at <http://www.medcomp.medicina.unipd.it/qparse/tool>

**Publication**: Michele Berselli, Enrico Lavezzo, Stefano Toppo, QPARSE: searching for long-looped or multimeric G-quadruplexes potentially distinctive and druggable, Bioinformatics, btz569, *10.1093/bioinformatics/btz569*

## **QPARSE_x.x.py**
### Requirements
QPARSE requires Python (version 2.7) and the numpy library.

To install numpy under unix environment (linux, osx):

	sudo pip install numpy

To get help with numpy installation check
[numpy docs](https://docs.scipy.org/doc/ "numpy documentation")

### Running QPARSE

To run QPARSE under unix environment (linux, osx):

	# move to QPARSE folder
	cd PATH/TO/QPARSE/FOLDER/

	# make QPARSE executable
	sudo chmod +x QPARSE_x.x.py

	# run QPARSE
	./QPARSE_x.x.py -i PATH/INPUT/FILE -o PATH/OUTPUT/FILE [OPTIONS]

	# !!! if the above is not working run
	python2.7 QPARSE_x.x.py -i PATH/INPUT/FILE -o PATH/OUTPUT/FILE [OPTIONS]


### Input
The program accepts as input both fasta and multi-fasta files.

### Parameters
The tool is very flexible and allows the user to select different parameters that can be defined to refine and customize the analysis. The parameters that can be used are:

***note**: defaults are shown in [], bool parameters do not require an argument. Parameters accepting numbers require integers.*

#### Basic parameters:
  - **-i**, **--inputsequence** PATH/TO/INPUTFILE: input file with sequence/s as fasta/multi-fasta.
  - **-o**, **--output** PATH/TO/OUTPUTFILE: output file to use for writing results.
  - **-b**, **--base** BASE [G]: letter that is used to detect the islands. G is used by default but other letters can be selected and used for the analysis instead (use UPPERCASE).
  - **-m**, **--minlen** MINLEN [2]: this parameter defines the length for the islands. It is also possible to select a range of lengths by adding the parameter **–M**, **--maxlen** MAXLEN to detect all the islands that are in the range of length [m..M]. While all the possible islands of the different lengths in the range are detected, only islands of the same length can be part of the same PQS. To allow for mismatched islands use the parameter **-x**, **-mismatch** MISMATCH [0] to set the maximum number of mismatched islands allowed.
  - **-L**, **--maxloop** MAXLOOP [12]: this parameter defines the maximum distance (loop distance) that is allowed between two consecutive islands within the same PQS.
  - **-g**, **--gapnum** GAPNUM [0]: this parameter defines the maximum number of gaps that can be opened per island. Alternatively, or in combination, the parameter **–l**, **--gaplen** GAPLEN [0] defines the maximum cumulative length of the gaps that is permitted per island. Finally, the parameter **–nocore** (bool) defines whether at least two consecutive bases are required to define an island. This Boolean parameter allows to detect also islands that do not have a “core” of at least two consecutive bases (e.g. G, GtGaG).
  - **-n**, **--islandnum** ISLANDNUM [4]: this parameter defines the number of islands that are required to be consecutive in the same PQS.
  - **-p**, **--perfect** PERFECT [1]: this parameter defines the minimum number of islands that are required to be “perfect” within each PQS. Since the tool allows to detect also islands that are degenerate, this parameter imposes a constraint for the minimum number of islands that cannot contain any degeneration within each pattern. If also PQS with no perfect islands are required use the parameter **-noperfect**, **--noperfect** (bool).
  -	**-fast**, **--fastermode** (bool): this parameter triggers a faster search mode that can be applied for islandnum > 4. Only regions with at least **--islandnum** islands are evaluated to build the graph, however, the graph is navigated searching for patterns of four islands to reduce the combinatorial complexity and speed up the analysis.

#### Loop symmetry check:
  -	**-sM**, **--simmetrymirror** (bool): evaluates the symmetry of the long loops (>= 7 nt) to improve the score, MIRROR symmetry is considered. Allows to detect longer loops with mirror properties that can form Hoogsteen-hairpins.
  -	**-sP**, **--simmetrypalindrome** (bool): evaluates the symmetry of the long loops (>= 7 nt) to improve the score, PALINDROMIC symmetry is considered. Allows to detect longer loops with palindromic properties that can form canonical-hairpins.
  -	**-sX**, **--simmetrymixed** (bool): evaluates the symmetry of the long loops (>= 7 nt) to improve the score, MIXED MIRROR-PALINDROMIC symmetry is considered. Allows to detect longer loops with mirror and palindromic properties that can form mixed-hairpins.
  -	**-sC**, **--simmetrycustom** PATH/TO/CUSTOM_MATRIX: evaluates the symmetry of the long loops (>= 7 nt) to improve the score, the custom substitution-matrix specified is used to test the alignment (see the template custom_substitution_matrix.txt).

#### Parameters to be used with caution:
  -	**-normal**, **--normalmode** (bool): searching for more than 12 islands **--fastermode** is used by default, this parameter override **--fastermode** and the graph is navigated searching for patterns of **--islandnum** islands *[caution!: the analysis can be computationally expensive]*.


### Examples of command lines
  - Perfect PQS with 4 islands of length 3 (-m), maximum loop distance of 7 nt (-L). Default search for G, -b C can be added to search for C instead:

		./QPARSE_x.x -i PATH/INPUT/FILE -o PATH/OUTPUT/FILE -m 3 -L 7 [-b C]

  - Perfect PQS with 4 islands of length in range [3..4] (-m, -M), maximum loop distance of 7 nt (-L). Default search for G, -b C can be added to search for C instead:

		./QPARSE_x.x -i PATH/INPUT/FILE -o PATH/OUTPUT/FILE -m 3 -M 4 -L 7 [-b C]

  - Degenerate PQS with 4 islands of length in range [3..4] (-m, -M), maximum loop distance of 7 nt (-L). Islands can have a maximum of 1 gap (-g) and a maximum gap length of 2 nt (-l). Default search for G, -b C can be added to search for C instead:

		./QPARSE_x.x -i PATH/INPUT/FILE -o PATH/OUTPUT/FILE -m 3 -M 4 -L 7 -g 1 -l 2  [-b C]

  - Degenerate PQS with 4 islands of length in range [3..4] (-m, -M), maximum loop distance of 7 nt (-L). Islands can have a maximum of 2 gap (-g) and a maximum gap length of 3 nt (-l). Also islands that do not have a “core” of at least two consecutive bases (-nocore) are detected (e.g. GcGtG). Default search for G, -b C can be added to search for C instead:

		./QPARSE_x.x -i PATH/INPUT/FILE -o PATH/OUTPUT/FILE -m 3 -M 4 -L 7 -g 2 -l 3 -nocore [-b C]

  - Degenerate PQS with 8 islands (-n) of length in range [3..4] (-m, -M), maximum loop distance of 7 nt (-L). Islands can have a maximum of 1 gap (-g) and a maximum gap length of 2 nt (-l). At least 5 non-degenerate islands (-p) are required. Default search for G, -b C can be added to search for C instead:

		./QPARSE_x.x -i PATH/INPUT/FILE -o PATH/OUTPUT/FILE -m 3 -M 4 -L 7 -g 1 -l 2 -n 8 -p 5 [-b C]

  - Degenerate PQS with 4 islands of length in range [3..4] (-m, -M), maximum loop distance of 15 nt (-L). Islands can have a maximum of 1 gap (-g) and a maximum gap length of 2 nt (-l). Long loops (>= 7 nt) are analyzed for mixed mirror-palindromic symmetrical properties (-sX). Default search for G, -b C can be added to search for C instead:

		./QPARSE_x.x -i PATH/INPUT/FILE -o PATH/OUTPUT/FILE -m 3 -M 4 -L 15 -g 1 -l 2 -sX [-b C]

### Output
The program returns at each position the highest scoring and more promising results.

#### Standard output

	#quadruplex     score   start   end     island_len
	>GeneID_1
	GGGG-tgt-GGGG-aca-GGGG-tgt-GGGG 40      3       27      4
	>GeneID_2
	GTGG--GGG--GGATG-taggt-GGG      16      3       22      3
	GTGG--GGG-ggat-GTAGG-t-GGG      16      3       22      3
	GTGG-g-GGG-gat-GTAGG-t-GGG      16      3       22      3
	GTGG-gg-GGG-at-GTAGG-t-GGG      16      3       22      3
	GGG--GGG-gat-GTAGG-t-GGG        19      5       22      3
	GGG-g-GGG-at-GTAGG-t-GGG        19      5       22      3
	GGG--GGG-at-GTAGG-t-GGG 19      6       22      3

The program returns a standard output that is structured in blocks. Each block corresponds to one of the sequences analyzed. Each block starts with a fasta-like header containing the information of the sequence as provided in the fasta input. The subsequent lines represent the detected PQS. Each line contains the sequence, the score, the start and the end indexes, and the length of the islands for the PQS. The values are separated by TAB. In the PQS sequence the islands are represented in CAPS LOCK and are connected to loops by '-': 'ISL-loop-ISL-loop-ISL-loop-...-ISL'. Null loops are represented as 'ISL--ISL'.

#### Output with loop symmetry check

	#quadruplex     score   start   end     island_len      loops_symmetry
	>GeneID_1
	GGCG-tta-GGG-aa-GGG-cgtcgaaagca-GGG     19      2       30      3       NA;NA;HmWuWm
	GGCG-tta-GGG-aa-GGG-cgtcgaaagcagggt-GGG 20      2       34      3       NA;NA;HlHHWHWll
	GGCG-tta-GGG-aagggcgtcgaaagca-GGG-t-GGG 19      2       34      3       NA;WWmmHWuHu;NA
	GGCG-ttagggaa-GGG-cgtcgaaagca-GGG-t-GGG 22      2       34      3       HmWW;HmWuWm;NA
	GGG-aa-GGG-cgtcgaaagca-GGG-t-GGG        22      9       34      3       NA;HmWuWm;NA
	>GeneID_2
	GGG-ac-GGG-ggccggc-GGG-ccac-GGG 27      3       27      3       NA;WWWu;NA
	GGG-acg-GGG-gccggc-GGG-ccac-GGG 30      3       27      3       NA;WWW;NA
	GGG-acgg-GGG-ccggc-GGG-ccac-GGG 24      3       27      3       NA;NA;NA

When the symmetry of the loops is considered, for each PQS an additional field is reported in the output. This field contains the information of the self-alignment for each of the loops. The self-alignment encodings for each loop are separated by ';'. In the encodings, 'W' identifies a Watson-Crick pairing, 'H' a Hoogsteen pairing, 'l-u' a gap opening and 'm' a mismatch. If the loop is shorter than 7 nt and the alignment is not calculated 'NA' is reported instead. This field is used to show the optimal alignment when using QPARSE_parser_x.x.py (see below).

## **QPARSE_parser_x.x.py**
Together with QPARSE, a Python script (Python v. 2.7) is also provided that can be used to better organize the raw output.

### Command line

	./QPARSE_parser_x.x.py -i PATH/INPUT/FILE -o PATH/OUTPUT/FILE [-g] [-m] [-x] [-s] [-a] [-maxLoop MAXLOOP]

### Parameters
***note**: bool parameters do not require an argument. Parameters accepting numbers require integers.*

#### General parameters:
  - **-i**, **--inputfile** PATH/TO/INPUTFILE: output file from QPARSE as input.
  - **-o**, **--outputfile** PATH/TO/OUTPUTFILE: file to store the formatted output.
  - **-g**, **--gff** (bool): by default the output of the parser is in tsv (TAB separated) format, this parameter converts the output to gff format.
  - **-m**, **--merge** (bool): using this parameter the overlapping PQS are merged to return the longest regions of consecutive islands whithin the loop distance.
  - **-x**, **--max** (bool): using this parameter only the maximum number of non overlapping PQS is returned.
  - **-maxLoop**, **--maxLoop** MAXLOOP: this parameter allows to define a maximum number of long loops (>= 7 nt) that is permitted in the results. The PQS with a higher number of long loops are filtered out.
  - **-s**, **--score** (bool): this parameter allows to order the PQS by score.

#### Parameters to use when loops symmetries are considered:
  - **-a**, **--alignment** (bool): this parameter allows to order the PQS by score, and shows the optimal alignment calculated for the linking loops.

#### Parameters to use with CAUTION (require mfold 3.6):
The software has been tested and is compatible with mfold version 3.6. Please refer to [mfold website](http://unafold.rna.albany.edu/?q=mfold "mfold resources") for resources and credits. mfold should be accessible from the command line (known to your PATH environment variable).

***WARNING**: mfold is computationally expensive, it is strongly recommended to use this utility only for a small number of sequences.*

  - **-mfold_s**, **--mfold_score** [DNA|RNA]: this parameter allows to order the PQS by score, in addition returns the energies calculated for the linking loops (>= 7 nt) using mfold. Specify DNA or RNA as parameter depending on your molecules.
  - **-mfold_a**, **--mfold_alignment** [DNA|RNA]: this parameter allows to order the PQS by score, and shows the more stable conformation and the energies calculated for the linking loops (>= 7 nt) using mfold. Specify DNA or RNA as parameter depending on your molecules.

### Output
#### Standard tsv (TAB separated) output

	#seqID  start   end     island_len      score   quadruplex
	GeneID_1     947     997     3       26      CCC-agccccctccgggccct-CCC-agcccctc-CCC-ttcctttccgcggc-CCC
	GeneID_2     913     944     3       26      CCC-cgccccgt-CCC-gacccct-CCC-gggtc-CCC

Each line contains the sequence ID as in the fasta input, the start and the end indexes, the length of the islands, the score and the sequence for the PQS. The values are separated by TAB. In the PQS sequence the islands are represented in CAPS LOCK and are connected to loops by '-': 'ISL-loop-ISL-loop-ISL-loop-...-ISL'. Null loops are represented as 'ISL--ISL'.

#### Output with alignments (-a)

	#seqID  start   end     island_len      score   quadruplex
	GeneID_1     947     997     3       26      CCC-agccccctccgggccct-CCC-agcccctc-CCC-ttcctttccgcggc-CCC
        Loop_1:
                agccccctc
                ||::||| :
                tcccggg-c
        Loop_2:
                -a-gc
                 | |:
                ctccc
        Loop_3:
                ttc-c-tt
                  | |  :
                cggcgcct
	GeneID_2     913     944     3       26      CCC-cgccccgt-CCC-gacccct-CCC-gggtc-CCC
        Loop_1:
                cgcc
                 :::
                tgcc
        Loop_2:
                gacc
                 |::
                -tcc
        Loop_3:
                --

                --

When showing the alignment, in between each PQS (in tsv format) it is shown the optimal alignment for each of the loops. '|' represent a Watson-Crick pairing, ':' represent a Hoogsteen pairing.

## **License**
Copyright (C) 2019 Michele Berselli

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

See http://www.gnu.org/licenses/ for more informations.

## **References**
**[1]** Characterization of G4–G4 Crosstalk in the c-KIT Promoter Region. Riccardo Rigo and Claudia Sissi. *Biochemistry* 2017.
**[2]** Formation of a Unique End-to-End Stacked Pair of G-Quadruplexes in the hTERT Core Promoter with Implications for Inhibition of Telomerase by G-Quadruplex-Interactive Ligands. SunMi L. Palumbo, Scot W. Ebbinghaus, and Laurence H. Hurley. *Journal of the American Chemical Society* 2009.
**[3]** Major G-Quadruplex Form of HIV-1 LTR Reveals a (3 + 1) Folding Topology Containing a Stem-Loop. Elena Butovskaya, Brahim Heddi, Blaž Bakalar, Sara N. Richter, and Anh Tuân Phan. *Journal of the American Chemical Society* 2018.
