# Snakemake workflow: Amplicon-SMF

## Authors
- Ben Doughty @bdoughty
- Jason Tan @yj-tan
- Julia Schaepe @juliaschaepe

## Overview:
The code here takes paired-end FASTQs from an amplicon SMF experiment and analyzes the single reads for methylation state. It returns both bulk measurements and single-molecule measurements, and the repo has additional code for downstream plotting. There are also scripts in `workflow/scripts` for calling binding (`classify_single_molecule_binding_v2.py`) and fitting the partition function model (`fit_partition_function_model_v3.py`).

## Requirements:
Snakemake handles the conda environments (two separate ones, for python2 and 3, see there for a specific list of packages and versions). So to run this, all you need is an environment containing Snakemake (7.15.2) and mamba (0.27.0). Memory requirements are usually pretty minimal (especially for data from a Miseq). Typically I run these on an interactive job on Stanford's HPC cluster with ~16G of memory, and it typically runs overnight.

## Input Files

### Samplesheet: 
The samplesheet contains all of the sample-level information you need to run. An example can be found in the directory below. Specifically, the columns it needs to include are:
- sample_name: English name for the sample, cannot be duplicated across the samplesheet (but you can make it up, and it doesn’t need to be the same as the FASTQ)
- fastq_R1: Full path to the R1 FASTQ file
- fastq_R2: Full path to the R2 FASTQ file
- amplicon_fa: Path to a FASTA file containing all of the amplicons you might want to align to, can be different for each sample. The sequences in the FASTAs can be longer than the reads, although then you will need to set “ignore_bounds” to TRUE
    - Note, these FASTAs should be the same strand as the read1 primer, i.e. if read1 primer binds “5’” of the sequence, then it should be the “top” strand on benchling (e.g. it should start with whatever is immediately downstream of R1).
    - However, because of the way we are designing primers, we only capture one actual strand, see below 
- experiment: English name for experiment-level grouping (i.e. all the samples with the same experiment will be output in a folder together
- filter_contigs: Necessary for backwards compatibility reasons, just say this is always TRUE (although this might change, as we update the read filtering code) [OPTIONAL; default TRUE]
- include_cpg: Whether or not to include CpGs in the analysis [OPTIONAL; default FALSE]
- no_endog_meth: Whether endogenous CpG methylation is expected to be prevalent or not (note that this is slightly different than the above, for instance, even if we didn’t include CpG MTase, there are some GpC sites that overlap CpG sites, and if we expect there to be a decent amount of endogenous methylation, we would have to ignore these) [OPTIONAL; default FALSE]
- ignore_bounds: Whether or not to mandate that all reads start/end within one “perfect-primer length” of the amplicons in amplicon_fa. This was included because sometimes there is slippage during PCR/sequencing that leads to reads that map (perfectly) in the middle of a molecule. The reason to turn this off would be to avoid having to duplicate all the amplicons for both E and P PCR in ExP, for example [OPTIONAL; default FALSE]
- read1_length: How long read1 (actually read2, since in pipeline read1 and read2 are swapped) was sequenced (note that I could have just gotten this from the FASTQs but I’m lazy, this will change) [OPTIONAL; default in config]
- read2_length: How long read2 (actually read1, since in pipeline read1 and read2 are swapped) was sequenced (note that I could have just gotten this from the FASTQs but I’m lazy, this will change) [OPTIONAL; default in config]
- bottom_strand: whether the strand captured by the SMF primers is the “bottom” strand (as it is for the TetO library) or not, tells this code which strand to search for conversion on [OPTIONAL; default TRUE]

### Configfile:
The configfile contains experiment-level information, specifically:
- samples: (Relative) path to samplesheet 
- alignment_score_fraction: how strict we should be on throwing away badly mapping reads (this and below are flags that affect how many reads we throw away vs. keep, when there are ambiguities, e.g. certain kinds of PCR errors or multimappers or chimeras, etc., the default I use is 0.8 which seems to work well across experiments)
- alignment_length_fraction: how strict we should be on throwing away reads that don’t map to enough of the library (this and above are flags that affect how many reads we throw away vs. keep, when there are ambiguities, e.g. certain kinds of PCR errors or multimappers or chimeras, etc., the default I use is 0.8 which seems to work well across experiments)
- read1_length: default read length for all samples (if you don’t want to specify for each)
- read2_length: default read length for all samples (if you don’t want to specify for each)
- unconverted_frac: fraction of non-GpC Cs that need to be converted for a read to count (0.85 seems to work well, except might be different transitioning into the genome)

Note: read1_length should actually be length of read2 from the sequencer since in the pipeline read1 and read2 are swapped. Similarly, read2_length should be the length of read1's from the sequencer.

## Running:
- First, pull a copy of the amplicon-smf code into the local directory:
`git clone https://github.com/GreenleafLab/amplicon-smf.git`
- Then all that’s left is to call Snakemake:
`snakemake -s amplicon-smf/workflow/Snakefile -c1 -k -w 15 --configfile [config.yaml]`
- Typically, the packages take 10s of minutes to install and the pipeline takes ~overnight to run.

## FAQs:
- Q: Why is bwameth failing to run for a subset of samples and telling me I need to index, didn’t you remember to include that step in the pipeline?
- A: As of Sept 27, 2022 this should be resolved. The issue previously was that bwameth was actually checking some auxillary index file when determining whether the fasta had been indexed already that we weren't explicitly requiring in the Snakemake rules. These auxillary files have since been added. Old answer from BD: This just happens, I have no idea why. Think it might have something to do with the modification times of various files being too close together? If you rerun the snakemake to figure out which failed and then manually reindex it seems to work. Sorry about that.
- Q: Why is everything "backwards"?
- A: This is kind of a subtle question having to do both with the intricacies of the bisulfite-aware aligners, but also because Ben is stupid. If you care, just reverse complement the FA you feed in :)

## Outputs:
Everything is output into the “results” folder. There is one subfolder per experiment. In each of the experiment folders, there is a folder called “plots”, which I look at first. It contains, per sample:
- A bulk_plots.pdf file, which has a bulk SMF trace and a methylation rate at GCs and other Cs (or alternately GCs, CGs, and other Cs if we are including CpGs) per amplicon in the FASTA
- A nonconverted_reads.png file, which shows the distribution of reads per amplicon in the sample, split by conversion status
- A nonconverted_reads.fraction.png file, which shows the distribution of conversion states (inferred from the non GpC Cs) per molecule and the cutoff used to filter reads (which I think is buried in that code but is somewhere \~85%)
- A wasted_reads.pdf file, which contains the info about how many reads were there to begin with and how many were retained. 
- A nuc_len_qc_plots.pdf, which tries to compute how long the longest streak of protection is per molecule, and fits some gaussian mixture models to figure out the average length of a nucleosome in the dataset. If there aren't that many TFs binding and messing things up, if the mode is around 130bp you're good!

Then, in each sample, there are basically only 2 things worth looking at:
- An amplicon_stats.txt file, which contains how many reads went to each library member and how many unique methylation patterns were observed, which gives a rough sense of the duplication rate
- In matrices/, there is one {sample}.{amplicon}.full_unclustered.matrix file, which contains a representation of all the reads (protected = 1, accessible = 0, and no info = -1). There’s also a png of a sample of 1,000 reads from that to visualize

There's also a QC folder, which currently only has fastqc results. I hope one day to collate some of the experimental stats per sample and output them all here. A work in progress.

## Downstream analyses:
### Running the binding model
`python workflow/scripts/classify_single_molecule_binding_v2.py --input /path/to/[samp]/matrices/[samp].[amp].full_unclustered.matrix --output /path/to/output/[samp].[amp].single_molecule_classification.txt --all_states_output /path/to/output/[samp].[amp].valid_states_table.txt --precomputed_states_file ./tmp/[amp].valid_states_table.pkl --plot /path/to/output/[samp].[amp].single_molecule_clustering.pdf --positions /path/to/positions.txt --amplicon_name [amp] --reads_to_plot 50`
The positions file can be created from the amplicon fasta file and a fasta file containing one entry per motif using `workflow/scripts/convert_fa_to_positions_for_script.py`.

### Running the partition function model
`python workflow/scripts/fit_partition_function_model_v3.py --basedir /path/to/output/ --sample [samp] --amplicons [comma,separated,list,of,amplicons] --plot /path/to/output/[samp].fit.pdf --model {'3param_nuc', '2param', '3param_tfcoop'} --output /path/to/output/[samp].fit.txt`

In amplicon-smf/workflow/scripts/ there are some other helpful downstream plotting files, although not a ton. As the project evolved the needs changed. I'm not gonna document them all here, but the titles are pretty informative. But basically they all just involved loading a matrix into pandas and doing plotting, so...

## To Do
- Reparallelize the bwameth step!
- Fix off by one error in matrix file generation (first base in reference is chopped off in matrix files). Likely problem is in convert_amplicon_fa_to_peaklist.py where start base is indicated as 1 but in python things are 0 indexed. So change this to a 0 and try, but this will cause downstream issues because of patchwork done in model, etc. around this indexing problem so wait until we can fix entire issue before touching this).
- Collate stats from the individual experiments and paste them in a table together.
