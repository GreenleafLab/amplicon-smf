# Snakemake workflow: Amplicon-SMF

## Authors
- Ben Doughty @bdoughty

## Overview:
The code here takes paired-end FASTQs from an amplicon SMF experiment and analyzes the single reads for methylation state. It returns both bulk measurements and single-molecule measurements, and the repo has additional code for downstream plotting.

## Requirements:
So far I have been running in my base conda environment, which has Snakemake, many scipy stack packages (e.g. numpy, scipy, matplotlib, seaborn, scikit-learn, etc.). I will try to make a Singularity container at some point, but it’s not ready yet. For now, you can try running it in your base environment or else you can make an environment based on the packages installed in mine (listed here: /oak/stanford/groups/wjg/bgrd/bin/211027_base_conda_env.txt). Then just get into an interactive job (sdev -p wjg -t 12:00:00 -m 32GB or so), navigate to the directory you want to run in, and follow the instructions below.

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
`git clone https://github.com/bdoughty/amplicon-smf.git`
- Then all that’s left is to call Snakemake:
`snakemake -s amplicon-smf/workflow/Snakefile -c1 -k -w 15 --configfile [config.yaml]`

Example directory: /oak/stanford/groups/wjg/bgrd/projects/smf/210817_smf_code_rewrite_reanalysis/211011_rerun_all (let me know if you can’t see this)

## FAQs:
- Q: Why is bwameth failing to run for a subset of samples and telling me I need to index, didn’t you remember to include that step in the pipeline?
- A: As of Sept 27, 2022 this should be resolved. The issue previously was that bwameth was actually checking some auxillary index file when determining whether the fasta had been indexed already that we weren't explicitly requiring in the Snakemake rules. These auxillary files have since been added. Old answer from BD: This just happens, I have no idea why. Think it might have something to do with the modification times of various files being too close together? If you rerun the snakemake to figure out which failed and then manually reindex it seems to work. Sorry about that.

## Outputs:
Everything is output into the “results” folder. There is one subfolder per experiment. In each of the experiment folders, there is a folder called “plots”, which I look at first. It contains, per sample:
- A bulk_plots.pdf file, which has a bulk SMF trace and a methylation rate at GCs and other Cs (or alternately GCs, CGs, and other Cs if we are including CpGs) per amplicon in the FASTA
- A nonconverted_reads.png file, which shows the distribution of reads per amplicon in the sample, split by conversion status
- A nonconverted_reads.fraction.png file, which shows the distribution of conversion states (inferred from the non GpC Cs) per molecule and the cutoff used to filter reads (which I think is buried in that code but is somewhere \~85%)
- A wasted_reads.pdf file, which contains the info about how many reads were included in downstream analyses (left) and how many filtered out (right). Note that the height of the left and right bars sums to the total number of reads (sorry, this is a bit confusing, I will change)

Then, in each sample, there are basically only 2 things worth looking at:
- An amplicon_stats.txt file, which contains how many reads went to each library member and how many unique methylation patterns were observed, which gives a rough sense of the duplication rate
- In matrices/, there is one {sample}.{amplicon}.full_unclustered.matrix file, which contains a representation of all the reads (protected = 1, accessible = 0, and no info = -1). There’s also a png of a sample of 1,000 reads from that to visualize

## Downstream analyses:
In amplicon-smf/workflow/scripts/ there are some other helpful downstream plotting files, although not a ton. You can check out /oak/stanford/groups/wjg/bgrd/projects/smf/210817_smf_code_rewrite_reanalysis/211011_rerun_all/further_analyses/log.sh for more information on how to run and the various outputs.

## To Do
- Make a Singularity container for these analyses (two, I guess, one for python2 and one for python3)
- Reparallelize the bwameth step!
- Fix off by one error in matrix file generation (first base in reference is chopped off in matrix files). Likely problem is in convert_amplicon_fa_to_peaklist.py where start base is indicated as 1 but in python things are 0 indexed. So change this to a 0 and try, but this will cause downstream issues because of patchwork done in model, etc. around this indexing problem so wait until we can fix entire issue before touching this).
