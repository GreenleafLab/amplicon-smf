from snakemake.utils import validate
import pandas as pd
from os import path
# from Bio import SeqIO


# # this container defines the underlying OS for each job when using the workflow
# # with --use-conda --use-singularity
# singularity: "docker://continuumio/miniconda3"

##### load config and sample sheets #####

# configfile: "amplicon-smf/config/config.yaml"
# validate(config, schema="../schemas/config.schema.yaml")

samplesheet = pd.read_csv(config["samples"], sep="\t").set_index("sample_name", drop=False)

# add back in the optional columns to the sample sheet, using defaults
if 'filter_contigs' not in samplesheet.columns:
    samplesheet['filter_contigs'] = True
if 'include_cpg' not in samplesheet.columns:
    samplesheet['include_cpg'] = False
if 'no_endog_meth' not in samplesheet.columns:
    samplesheet['no_endog_meth'] = False
if 'ignore_bounds' not in samplesheet.columns:
    samplesheet['ignore_bounds'] = False
if 'read1_length' not in samplesheet.columns:
    samplesheet['read1_length'] = config['read1_length']
if 'read2_length' not in samplesheet.columns:
    samplesheet['read2_length'] = config['read2_length']
if 'bottom_strand' not in samplesheet.columns:
    samplesheet['bottom_strand'] = True
if 'deaminase' not in samplesheet.columns:
    samplesheet['deaminase'] = False

# validate(samples, schema="../schemas/samples.schema.yaml")

def all_input(wildcards):
    wanted_input = []

    for sample in samplesheet['sample_name']:
        # amplicon_fa = samplesheet.loc[sample, 'amplicon_fa']
        experiment = samplesheet.loc[sample, 'experiment']

        wanted_input.append('results/{e}/plots/{s}.bulk_plots.pdf'.format(e=experiment, s=sample))
        wanted_input.append('results/{e}/plots/{s}.bulk_plots_from_matrices.pdf'.format(e=experiment, s=sample))
        if True: #samplesheet.loc[sample, 'filter_contigs']:
            wanted_input.append('results/{e}/plots/{s}.wasted_reads.pdf'.format(e=experiment, s=sample))
        wanted_input.append('results/{e}/{s}/{s}.amplicon_stats.txt'.format(e=experiment, s=sample))
        wanted_input.append('results/{e}/plots/{s}.nuc_len_qc_plots.pdf'.format(e=experiment, s=sample))
        wanted_input.append('results/qc/fastqc/fastqc.txt')
        # wanted_input.append('results/qc/fastqc/{}_fastqc.html'.format(path.basename(samplesheet.loc[sample,'fastq_R1']).replace('.fastq.gz','')))
        # wanted_input.append('results/qc/fastqc/{}_fastqc.html'.format(path.basename(samplesheet.loc[sample,'fastq_R2']).replace('.fastq.gz','')))

    # wanted_input.append('results/qc/fastqc/fastqc.txt')

    return wanted_input
