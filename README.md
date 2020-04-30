# chromeister
An ultra fast, heuristic approach to detect conserved signals in extremely large pairwise genome comparisons.

## Requirements

GCC compiler (any version that is not completely outdated should do) and the R programming language with the package "dplyr".

Simply download the .zip and unzip it, or clone the repository.
Then issue the following command:

1. cd chromeister/src && make all

And then open an R session and install the dplyr package by doing:

1. R
2. install.packages("dplyr")

This should install the R package dplyr.

If the installation finished without errors, you are ready to go! In case you can not install dplyr (or do not wish to) you can still use CHROMEISTER (but the dotplot will show no grid).

## Use

There are several ways in which CHROMEISTER can be used. The simplest one is to run a 1-vs-1 comparison and then compute the score and the plot.
To do so, use the binaries at the bin folder:

### Simple execution

You can run CHROMEISTER directly by issuing:

CHROMEISTER -query seqX -db seqY -out dotplot.mat && Rscript compute_score.R dotplot.mat 1000

If you do not want a grid on the output dotplot (which is recommended when running comparisons with a lot of scaffolds for instance) then run the same command but replace compute_score by compute_score-nogrid, see below:

CHROMEISTER -query seqX -db seqY -out dotplot.mat && Rscript compute_score-nogrid.R dotplot.mat 1000

The 1000 value is the default size of dotplot.mat, i.e. the resolution of the matrix -- if you want to change this (for example to generate a larger image (if you use 2000 it will generate a plot of 2000x2000, so be careful) include also the parameter -dimension in CHROMEISTER. Example command with larger resolution:

CHROMEISTER -query seqX -db seqY -out dotplot.mat -dimension 2000 && Rscript compute_score.R dotplot.mat 2000

or you can also use the script that is in the bin folder (which will do the above for you):

run_and_plot_chromeister.sh (input sequence A) (input sequence B) (KMER size) (DIMENSION of plot) (inexactitude level) [optional: grid]

(see parameters at the end) (the grid keyword at the end can be included/omitted depending if you want grid in the output dotplot)

This will generate the following items:

*	Comparison matrix, i.e. a scaled matrix containing the unique and inexact hits
*	Plot of the comparison with the automatic scoring distance and grid separating different sequences (chromosomes for instance)
*	CSV file containing the coordinates of each sequence/chromosome contained within the query and the reference
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


## Converting CHROMEISTER signal into alignments

First of all, consider whether it is interesting or not to use CHROMEISTER for "fine-grained" results. CHROMEISTER is recommended for VERY coarse-grained and full-genome comparisons in order to quickly assess similarity between genomes. Thus it does NOT produce alignments. However, if you find yourself in a situation where you want to convert the signal of CHROMEISTER into alignments (e.g. two large genomes), this can be done. The following tutorial shows how to do it, with human chromosome X and mouse chromosome X as example:

1. First, run CHROMEISTER like this:

	./CHROMEISTER -query HOMSA.Chr.X.fasta -db MUSMU.Chr.X.fasta -out dotplot.mat -dimension 1000 && Rscript compute_score.R dotplot.mat 1000

2. Check the "dotplot.mat.filt.png" corresponding to the dotplot between both chromosomes to see if there is any similarity. If so, proceed to next step.

3. Clone the following repository: https://github.com/estebanpw/gecko

	git clone https://github.com/estebanpw/gecko

4. Switch branch to the one named "inmemory_guided_chrom" and compile it. To do so, issue the following command:
	
	cd gecko && git checkout inmemory_guided_chrom && make all -C src
	

5. Now run the script "guidefastas" in the bin folder. See below:

	bin/guidefastas.sh HOMSA.Chr.X.fasta MUSMU.Chr.X.fasta hits-XY-dotplot.mat.hits 1000 100 60 32

Note (1): remember to include the full path to the sequences.
Note (2): the "hits-XY-dotplot.mat.hits" file is produced by CHROMEISTER in step 1. Copy it to the folder or include full path.
Note (3): the parameters following in the command "1000 200 75 32" are namely (1) size of dotplot, (2) minimum length that an alignment must have to be reported, (3) minimum similarity from 0-100, (4) k-mer seed size (use 32 for chromomsome-like sequences).

This step can take several minutes, e.g. using 1 CPU this execution took around 9-10 minutes.

6. A CSV file containing the alignments coordinates can be found in the folder all-results/master.csv. You can download it here if you wish to do so: http://mango.ac.uma.es/compartir/HOMSA_X-MUSMU_X.csv

7. If you also wish to visually contrast annotations to the alignments, you can use our genomic browser at https://pistacho.ac.uma.es/. To do so just follow the user guide available at https://pistacho.ac.uma.es/static/data/GeckoMGV-UserGuide.pdf


## Parameters

USAGE:
* -query: 	sequence A in fasta format
* -db: 		sequence B in fasta format
* -out:		output matrix
* -kmer       Integer:   k>1 (default 32) Use 32 for chromosomes and genomes and 16 for small bacteria
* -diffuse    Integer:   z>0 (default 4) Use 4 for everything - if using large plant genomes you can try using 1
* -dimension  Size of the output matrix and plot. Integer:   d>0 (default 1000) Use 1000 for everything that is not full genome size, where 2000 is recommended













