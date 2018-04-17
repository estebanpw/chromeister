# chromeister
An ultra fast, heuristic approach to detect conserved signals in extremely large pairwise genome comparisons.

## Requirements

GCC compiler (any version that is not completely outdated should do) and the R-base package (default R installation should do) with no extra packages.
Simply download the .zip and unzip it, or clone the repository (currently the active branch is "histo-kmer-ultra" DO NOT USE ANY OTHER).
Then issue the following commands:

cd chromeister/src && make all

You are ready to go!

## Use

There are several ways in which CHROMEISTER can be used. The simplest one is to run a 1-vs-1 comparison and then compute the score and the plot.
To do so, use the binaries at the bin folder:

CHROMEISTER -query seqX -db seqY -out dotplot.mat && Rscript compute_score.R dotplot.mat

This will generate the comparison matrix, the plot of the comparison with the automatic scoring and the guides to be used in an exhaustive GECKO comparison.

More to be added soon!













