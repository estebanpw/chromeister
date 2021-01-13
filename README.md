[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/chromeister/README.html)
[![European Galaxy server](https://img.shields.io/badge/usegalaxy-.eu-brightgreen?logo=data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABgAAAASCAYAAABB7B6eAAAABGdBTUEAALGPC/xhBQAAACBjSFJNAAB6JgAAgIQAAPoAAACA6AAAdTAAAOpgAAA6mAAAF3CculE8AAAACXBIWXMAAAsTAAALEwEAmpwYAAACC2lUWHRYTUw6Y29tLmFkb2JlLnhtcAAAAAAAPHg6eG1wbWV0YSB4bWxuczp4PSJhZG9iZTpuczptZXRhLyIgeDp4bXB0az0iWE1QIENvcmUgNS40LjAiPgogICA8cmRmOlJERiB4bWxuczpyZGY9Imh0dHA6Ly93d3cudzMub3JnLzE5OTkvMDIvMjItcmRmLXN5bnRheC1ucyMiPgogICAgICA8cmRmOkRlc2NyaXB0aW9uIHJkZjphYm91dD0iIgogICAgICAgICAgICB4bWxuczp0aWZmPSJodHRwOi8vbnMuYWRvYmUuY29tL3RpZmYvMS4wLyI+CiAgICAgICAgIDx0aWZmOlJlc29sdXRpb25Vbml0PjI8L3RpZmY6UmVzb2x1dGlvblVuaXQ+CiAgICAgICAgIDx0aWZmOkNvbXByZXNzaW9uPjE8L3RpZmY6Q29tcHJlc3Npb24+CiAgICAgICAgIDx0aWZmOk9yaWVudGF0aW9uPjE8L3RpZmY6T3JpZW50YXRpb24+CiAgICAgICAgIDx0aWZmOlBob3RvbWV0cmljSW50ZXJwcmV0YXRpb24+MjwvdGlmZjpQaG90b21ldHJpY0ludGVycHJldGF0aW9uPgogICAgICA8L3JkZjpEZXNjcmlwdGlvbj4KICAgPC9yZGY6UkRGPgo8L3g6eG1wbWV0YT4KD0UqkwAAAn9JREFUOBGlVEuLE0EQruqZiftwDz4QYT1IYM8eFkHFw/4HYX+GB3/B4l/YP+CP8OBNTwpCwFMQXAQPKtnsg5nJZpKdni6/6kzHvAYDFtRUT71f3UwAEbkLch9ogQxcBwRKMfAnM1/CBwgrbxkgPAYqlBOy1jfovlaPsEiWPROZmqmZKKzOYCJb/AbdYLso9/9B6GppBRqCrjSYYaquZq20EUKAzVpjo1FzWRDVrNay6C/HDxT92wXrAVCH3ASqq5VqEtv1WZ13Mdwf8LFyyKECNbgHHAObWhScf4Wnj9CbQpPzWYU3UFoX3qkhlG8AY2BTQt5/EA7qaEPQsgGLWied0A8VKrHAsCC1eJ6EFoUd1v6GoPOaRAtDPViUr/wPzkIFV9AaAZGtYB568VyJfijV+ZBzlVZJ3W7XHB2RESGe4opXIGzRTdjcAupOK09RA6kzr1NTrTj7V1ugM4VgPGWEw+e39CxO6JUw5XhhKihmaDacU2GiR0Ohcc4cZ+Kq3AjlEnEeRSazLs6/9b/kh4eTC+hngE3QQD7Yyclxsrf3cpxsPXn+cFdenF9aqlBXMXaDiEyfyfawBz2RqC/O9WF1ysacOpytlUSoqNrtfbS642+4D4CS9V3xb4u8P/ACI4O810efRu6KsC0QnjHJGaq4IOGUjWTo/YDZDB3xSIxcGyNlWcTucb4T3in/3IaueNrZyX0lGOrWndstOr+w21UlVFokILjJLFhPukbVY8OmwNQ3nZgNJNmKDccusSb4UIe+gtkI+9/bSLJDjqn763f5CQ5TLApmICkqwR0QnUPKZFIUnoozWcQuRbC0Km02knj0tPYx63furGs3x/iPnz83zJDVNtdP3QAAAABJRU5ErkJggg==)](https://usegalaxy.eu/root?tool_id=chromeister)

# chromeister
An ultra fast, heuristic approach to detect conserved signals in extremely large pairwise genome comparisons.

## Requirements

GCC compiler (any version that is not completely outdated should do), the R programming language (base installation) and python3 (tested on 3.5 and 3.8). Please make sure that on your linux CC resolves to GCC, otherwise it might not work.

Simply download the .zip and unzip it, or clone the repository.
Then issue the following command:

`cd chromeister && make all -C src/ && python3 -m venv chromeisterenv && source chromeisterenv/bin/activate && pip install -r src/requirements.txt`

This will compile CHROMEISTER and create a virtualenv where the python libraries will be installed (see `src/requirements.txt`)

If the installation finished without errors, you are ready to go! If you encounter the following problem `ImportError: No module named 'skbuild'` you might need to do `source chromeisterenv/bin/activate && pip install --upgrade pip` and then run `pip install -r src/requirements.txt` to finish the installation. 

**NOTE**: python and its libraries are only used for the detection of events. If the binaries compile (i.e. the `make all`) then you can still run CHROMEISTER and plot the results, even if the python installation did not work.

## Use

There are several ways in which CHROMEISTER can be used. The simplest one is to run a 1-vs-1 comparison and then compute the score and the plot.
To do so, use the binaries at the bin folder:

### Simple execution

You can run CHROMEISTER directly by issuing:

```CHROMEISTER -query seqX -db seqY -out dotplot.mat && Rscript compute_score.R dotplot.mat 1000```

If you do not want a grid on the output dotplot (which is recommended when running comparisons with a lot of scaffolds for instance) then run the same command but replace compute_score by compute_score-nogrid, see below:

```CHROMEISTER -query seqX -db seqY -out dotplot.mat && Rscript compute_score-nogrid.R dotplot.mat 1000```

The 1000 value is the default size of dotplot.mat, i.e. the resolution of the matrix -- if you want to change this (for example to generate a larger image (if you use 2000 it will generate a plot of 2000x2000, so be careful) include also the parameter -dimension in CHROMEISTER. Example command with larger resolution:

```CHROMEISTER -query seqX -db seqY -out dotplot.mat -dimension 2000 && Rscript compute_score.R dotplot.mat 2000```

And if you want to run the events detection, use (make sure that your virtualenv `chromeisterenv` is in the chromeister root folder:

`source ../chromeisterenv/bin/activate && python3 bin/detect_events.py dotplot.mat.raw.txt`

This will generate a `dotplot.mat.events.txt` file containing the detected events and classified. If you want to get a plot of the signal with the overlapped detected events, issue the same command but add at the end the parameter `png` (separated by a space).

You can also use the script that is in the bin folder (which will do all of the above for you):

```run_and_plot_chromeister.sh (input sequence A) (input sequence B) (KMER size) (DIMENSION of plot) (inexactitude level) [optional: grid]```

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
	 ```allVsAll.sh <sequences folder> <extension (e.g. fasta)> <matrix size (1000 for chromosomes, 2000 for full genomes)> <kmer size 32 (32 is best)> <inexactitude level (4 is recommended)> ```

* Comparing two folders containing sequences. This accounts for n * m comparisons, therefore it will compare ALL to ALL. Use this for instance to compare all chromosomes of one genome to all chromosomes of another genome.
	* To run this mode, use the script in the bin folder:
	 ```allVsAll_incremental.sh <sequences folder 1> <sequences folder 2> <extension (e.g. fasta)> <matrix size (1000 for chromosomes, 2000 for full genomes)> <kmer size 32 (32 is best)> <inexactitude level (4 is recommended)>```

At the end of both comparisons, an index will be created summarizing the scores per each comparison. This index has the following format (see header and example below):
header: <SpX, SpY, IDX, IDY, IMG, CHNumberX, CHNumberY, Score, LengthX, LengthY>
example: BRAOL.Chr.C1,BRAOL.Chr.C2,>C1 dna:chromosome chromosome:v2.1:C1:1:43764888:1 REF,>C2 dna:chromosome chromosome:v2.1:C2:1:52886895:1 REF,BRAOL.C
hr.C1.fasta-BRAOL.Chr.C2.fasta.mat.filt.png,C1,C2, 0.996,43764888,52886895


Notice that you can easily run this in parallel by just re-issuing the command (i.e. execute same command as many times as you want, each time another core will help in the processing).


## Converting CHROMEISTER signal into alignments

First of all, consider whether it is interesting or not to use CHROMEISTER for "fine-grained" results. CHROMEISTER is recommended for VERY coarse-grained and full-genome comparisons in order to quickly assess similarity between genomes. Thus it does NOT produce alignments. However, if you find yourself in a situation where you want to convert the signal of CHROMEISTER into alignments (e.g. two large genomes), this can be done. The following tutorial shows how to do it, with human chromosome X and mouse chromosome X as example:

1. First, run CHROMEISTER like this:

	```./CHROMEISTER -query HOMSA.Chr.X.fasta -db MUSMU.Chr.X.fasta -out dotplot.mat -dimension 1000 && Rscript compute_score.R dotplot.mat 1000```

2. Check the "dotplot.mat.filt.png" corresponding to the dotplot between both chromosomes to see if there is any similarity. If so, proceed to next step.

3. Clone the following repository: https://github.com/estebanpw/gecko

	```git clone https://github.com/estebanpw/gecko```

4. Switch branch to the one named "inmemory_guided_chrom" and compile it. To do so, issue the following command:
	
	```cd gecko && git checkout inmemory_guided_chrom && make all -C src```
	

5. Now run the script "guidefastas" in the bin folder. See below:

	```bin/guidefastas.sh HOMSA.Chr.X.fasta MUSMU.Chr.X.fasta hits-XY-dotplot.mat.hits 1000 100 60 32```
	
	You can add the keyword `alignments` and/or `names` at the end of the command. The first one will extract the alignments and show them on the screen (you can redirect them with `>`). The second one will output the names of the sequences to which each fragment belongs (as opposed to numbers). If you want to show the alignments and have names of sequences as well, run it as follows:

	```bin/guidefastas.sh HOMSA.Chr.X.fasta MUSMU.Chr.X.fasta hits-XY-dotplot.mat.hits 1000 100 60 32 alignments names```

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

## Test data

You can test CHROMEISTER with the two mycoplasma sequences provided in the 'test-data' folder. You can do so by running the following commands (from within the test-data folder):

```../bin/CHROMEISTER -query mycoplasma-232.fasta -db mycoplasma-7422.fasta -out mycoplasma-232-7422.mat -dimension 500```
```Rscript ../bin/compute_score.R mycoplasma-232-7422.mat 500```

Note: in this example we used size 500 since the two sequences are quite small.

## Example runs

### Chromosome example 

Comparing two chromosomes (Homo sapiens Chr X vs Mus musculus Chr X) in a minute:

`run_and_plot_chromeister.sh HOMSA.Chr.X.fasta MUSMU.Chr.X.fasta 32 1000 4`

![HOMSA vs MUSMU](https://github.com/estebanpw/chromeister/blob/master/images/HOMSA.Chr.X.fasta-MUSMU.Chr.X.fasta.mat.filt.png)


### Multi-fasta

Comparing the full genome of the Gallus gallus against the Meleagris gallopavo with all their chromosomes including auto-generated grid:

`run_and_plot_chromeister.sh GALGA.Chr.complete.fasta MELGA.Chr.complete.fasta 32 1000 4 grid`

![MULTI-FASTA](https://github.com/estebanpw/chromeister/blob/master/images/GALGA.Chr.complete.fasta-MELGA.Chr.complete.fasta.mat.filt.png)

### Events detection

From the `Homo sapiens Chr X vs Mus musculus Chr X` comparison you can inspect the file `HOMSA.Chr.X.fasta-MUSMU.Chr.X.fasta.mat.events.txt` which contains the classified events. You can also plot the overlap between detected events and original plot by doing:

`python3 chromeister/bin/detect_events.py HOMSA.Chr.X.fasta-MUSMU.Chr.X.fasta.mat.raw.txt png`

which will generate the following plot:

![DETECTED-EVENTS](https://github.com/estebanpw/chromeister/blob/master/images/HOMSA.Chr.X.fasta-MUSMU.Chr.X.fasta.mat.events.png)

Notice that the orange color represents detected events whereas green represents undetected blocks. The detection is based on the Houghs transform and therefore there are some parameters that can be used to change the minimum length to detect an event (e.g. if we want single points as well) or the minimum gap to join two blocks. These parameters can be changed in the `detect_events.py` file. The current configuration is mostly tailored to detect larger events and no single points.

### Fine-grained run

Comparing some mycoplasma hyopneumoniae genomes (not all has to be mammalian or plant chromosomes!).

![SMALL-BACTERIA](https://github.com/estebanpw/chromeister/blob/master/images/NC_014.fasta-NC_017.fasta.mat.filt.png)

This is afterwards ran with the `inmemory_guided_chrom` branch of GECKO (see [here](https://github.com/estebanpw/gecko/tree/inmemory_guided_chrom)) using the following command line (supply the CHROMEISTER hits file):

`gecko/bin/guidefastas.sh NC_014.fasta NC_017.fasta hits-XY-NC_014.fasta-NC_017.fasta.mat.hits 1000 50 60 32 names`

and we get a csv `NC_014-NC_017.csv` which includes:

```
Frag,1218913,3918445,1221077,3916281,r,0,2165,8644,2163,99.82,1.00,_gi|304372805|ref|NC_014448.1|,_gi|385858114|ref|NC_017519.1|
[...]
Frag,6618091,5316746,6620403,5319058,f,0,2313,8524,2222,92.13,0.96,_gi|321309518|ref|NC_014970.1|,_gi|385858893|ref|NC_017520.1|
[...]
Frag,3062371,6264699,3064642,6266970,f,0,2272,9080,2271,99.91,1.00,_gi|313664890|ref|NC_014751.1|,_gi|392388518|ref|NC_017521.1|
```
Where the last two columns indicate the origin of the fragment. For instance, in this case, the first fragment belongs to `gi|304372805|ref|NC_014448.1|` and `gi|385858114|ref|NC_017519.1|` which are sequences `2` and `5` in the `x` and `y` axis in the previous plot, respectively.

Remember that your `csv` file can be uploaded and interactively inspected [here](https://pistacho.ac.uma.es/)!

## Help

1. Hanging output and program does not finish
	If you experience this kind of output:

	[INFO] Generating a 1000x1000 matrix

	[INFO] Loading database

	100%...[INFO] Database loaded and of length 70039485.

	[INFO] Ratios: Q [6.658880e+04] D [7.003949e+04]. Lenghts: Q [66588797] D [70039485]

	[INFO] Pixel size: Q [1.501754e-05] D [1.427766e-05].

	[INFO] Computing absolute hit numbers.

	100%...Scanning hits table.

	100%...

	[INFO] Query length 66588797.

	[INFO] Writing matrix.

	[INFO] Found 238693 unique hits for z = 4.

	But the program doesnt finish (it "hangs"), then your input sequences probably contain a lot of sequences (i.e. a multifasta with hundreds of contigs). To fix this, simply run the same command but instead of using the "compute_score.R" use the "compute_score-nogrid.R" script. This will remove the drawing of the grid which can get overflown when using too many sequences.

## Citing

If you use or have used CHROMEISTER in your research, please cite the following article:

Esteban Pérez-Wohlfeil, Sergio Diaz-del-Pino and Oswaldo Trelles. "Ultra-fast genome comparison for large-scale genomic experiments." Scientific reports 9, no. 1 (2019): 1-10.

[Link to manuscript](https://www.nature.com/articles/s41598-019-46773-w)




