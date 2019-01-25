# chromeister
An ultra fast, heuristic approach to detect conserved signals in extremely large pairwise genome comparisons.

## Requirements

GCC compiler (any version that is not completely outdated should do) and the R-base package (default R installation should do) with no extra packages.
Simply download the .zip and unzip it, or clone the repository.
Then issue the following commands:

cd chromeister/src && make all

You are ready to go!

## Use

There are several ways in which CHROMEISTER can be used. The simplest one is to run a 1-vs-1 comparison and then compute the score and the plot.
To do so, use the binaries at the bin folder:

### Simple execution

You can run CHROMEISTER directly by issuing:

CHROMEISTER -query seqX -db seqY -out dotplot.mat && Rscript compute_score.R dotplot.mat

or you can also use the script that is in the bin folder (which will do the above for you):

run_and_plot_chromeister.sh (input sequence A) (input sequence B) (KMER size) (DIMENSION of plot) (inexactitude level)

(see parameters at the end)

This will generate the following items:

*	Comparison matrix, i.e. a scaled matrix containing the unique and inexact hits
*	Plot of the comparison with the automatic scoring distance
*	Events file. A text file where each row is a synteny block. Note: these events are Large-Scale Genome Rearrangements heuristically determined and classified as {Synteny block, transposition, inversion, ...} - but this is only an informative labelling that only considers coordinates - do not blindly believe in the classification, but rather do your own labelling based on the events.
*	Guides to be used in an exhaustive GECKO comparison (reduces runtime)

### All vs All execution

You can run massive all versus all comparisons in two diferent ways:

* Comparing all the sequences in one folder. This accounts for 1/2 * n * (n+1) comparisons, hence it will not compare sequence B to sequence A if the comparison for sequence A to sequence B already existed.
	* To run this mode, use the script in the bin folder:
	 allVsAll.sh <sequences folder> <extension (e.g. fasta)> <matrix size (1000 for chromosomes, 2000 for full genomes)> <kmer size 32 (32 is best)> <inexactitude level (4 is recommended)> 

* Comparing two folders containing sequences. This accounts for n * m comparisons, therefore it will compare ALL to ALL. Use this for instance to compare all chromosomes of one genome to all chromosomes of another genome.
	* To run this mode, use the script in the bin folder:
	 allVsAll_incremental.sh <sequences folder 1> <sequences folder 2> <extension (e.g. fasta)> <matrix size (1000 for chromosomes, 2000 for full genomes)> <kmer size 32 (32 is best)> <inexactitude level (4 is recommended)>

At the end of both comparisons, an index will be created summarizing the scores per each comparison. This index has the following format (see header and example below):
header: <SpX, SpY, IDX, IDY, IMG, CHNumberX, CHNumberY, Score, LengthX, LengthY>
example: BRAOL.Chr.C1,BRAOL.Chr.C2,>C1 dna:chromosome chromosome:v2.1:C1:1:43764888:1 REF,>C2 dna:chromosome chromosome:v2.1:C2:1:52886895:1 REF,BRAOL.C
hr.C1.fasta-BRAOL.Chr.C2.fasta.mat.filt.png,C1,C2, 0.996,43764888,52886895


Notice that you can easily run this in parallel by just re-issuing the command (i.e. execute same command as many times as you want, each time another core will help in the processing).



## Parameters

USAGE:
* -query: 	sequence A in fasta format
* -db: 		sequence B in fasta format
* -out:		output matrix
* -kmer       Integer:   k>1 (default 32) Use 32 for chromosomes and genomes and 16 for small bacteria
* -diffuse    Integer:   z>0 (default 4) Use 4 for everything - if using large plant genomes you can try using 1
* -dimension  Size of the output matrix and plot. Integer:   d>0 (default 1000) Use 1000 for everything that is not full genome size, where 2000 is recommended













