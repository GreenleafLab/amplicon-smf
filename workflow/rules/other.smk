rule fastqc:
    input:
        lambda wildcards: samplesheet[['fastq_R1','fastq_R2']].to_numpy().flatten().tolist()
    output:
        'results/qc/fastqc/fastqc.txt'
    conda:
        "envs/python3_v7.yaml"
    params:
        outdir='results/qc/fastqc'
    shell:
        'fastqc -o {params.outdir} -f fastq {input}; touch {output}'

# rule fastqc:
#     input:
#         ''
#     output:
#         'results/qc/fastqc/{read}_fastqc.html'
#     conda:
#         "envs/python3_v7.yaml"
#     params:
#         outdir='results/qc/fastqc'
#     shell:
#         'fastqc -o {params.outdir} -f fastq {input}; touch {output}'

def get_fastq(wildcards):
    # bottom_strand = samplesheet.loc[wildcards.sample,'bottom_strand']
    # if bottom_strand:
    #     fastq_col = 'fastq_R{}'.format(3-int(wildcards.read))
    # else:
    #     fastq_col = 'fastq_R{}'.format(wildcards.read)
    # fastq_col = 'fastq_R{}'.format(wildcards.read)
    fastq_col = 'fastq_R{}'.format(3-int(wildcards.read))
    fastq = samplesheet.loc[wildcards.sample, fastq_col]
    return fastq

rule reverse_complement_fastq:
    input:
        fq=get_fastq#,
        #fa='results/{experiment}/{sample}/tmp/{sample}.amplicon.revcomp.fa'
    output:
        # 'results/{sample}/tmp/{sample}.read{read}.revcomp.fastq.gz'
        temp('results/{experiment}/{sample}/tmp/{sample}_revcomp_R{read}_001.fastq.gz')
    conda:
        "envs/python3_v7.yaml"
    shell:
        # 'ml fastx_toolkit/0.0.14; zcat {input.fq} | fastx_reverse_complement -i - -Q33 -z > {output}'
        'cp {input.fq} {output}'

def get_amplicon(wildcards):
    return samplesheet.loc[wildcards.sample, "amplicon_fa"]

rule reverse_complement_fasta:
    input:
        get_amplicon
    output:
        'results/{experiment}/{sample}/tmp/{sample}.amplicon.revcomp.fa'
    params:
        rc_or_not=lambda wildcards: 'fastx_reverse_complement -i - -o' if samplesheet.loc[wildcards.sample,'bottom_strand'] else 'cat > '
    conda:
        "envs/python3_v7.yaml"
    shell:
        # "ml fastx_toolkit/0.0.14; awk '/^>/ {{print($0)}}; /^[^>]/ {{print(toupper($0))}}' {input} | fastx_reverse_complement -i - -o {output}"
        # "ml fastx_toolkit/0.0.14; awk '/^>/ {{print($0)}}; /^[^>]/ {{print(toupper($0))}}' {input} | {params.rc_or_not} {output}"
        "awk '/^>/ {{print($0)}}; /^[^>]/ {{print(toupper($0))}}' {input} | {params.rc_or_not} {output}"

# rule reverse_complement_fasta:
#     input:
#         get_amplicon
#     output:
#         "results/{experiment}/{sample}/tmp/{sample}.amplicon.revcomp.fa"
#     params:
#         # 1 = reverse-complement, 0 = just uppercase
#         rc_flag=lambda wildcards: 1 if samplesheet.loc[wildcards.sample, "bottom_strand"] else 0
#     conda:
#         "envs/python3_v7.yaml"
#     shell:
#         r"""
#         awk -v rc={params.rc_flag} '
#             BEGIN {{
#                 comp["A"]="T"; comp["T"]="A"; comp["G"]="C"; comp["C"]="G";
#                 comp["a"]="T"; comp["t"]="A"; comp["g"]="C"; comp["c"]="G";
#                 comp["N"]="N"; comp["n"]="N";
#             }}
#             /^>/ {{
#                 # If we already collected a previous record, flush it
#                 if (NR > 1) {{
#                     if (rc == 1) {{
#                         # reverse complement, uppercase
#                         n = length(seq);
#                         rev = "";
#                         for (i = n; i >= 1; i--) {{
#                             b = toupper(substr(seq, i, 1));
#                             if (b in comp) rev = rev comp[b];
#                             else rev = rev "N";
#                         }}
#                         print header;
#                         print rev;
#                     }} else {{
#                         # just uppercase
#                         print header;
#                         print toupper(seq);
#                     }}
#                 }}
#                 header = $0;
#                 seq = "";
#                 next;
#             }}
#             /^[^>]/ {{
#                 seq = seq $0;
#                 next;
#             }}
#             END {{
#                 if (header != "") {{
#                     if (rc == 1) {{
#                         n = length(seq);
#                         rev = "";
#                         for (i = n; i >= 1; i--) {{
#                             b = toupper(substr(seq, i, 1));
#                             if (b in comp) rev = rev comp[b];
#                             else rev = rev "N";
#                         }}
#                         print header;
#                         print rev;
#                     }} else {{
#                         print header;
#                         print toupper(seq);
#                     }}
#                 }}
#             }}
#         ' {input} > {output}
#         """

rule index_fasta:
    input:
        fa='results/{experiment}/{sample}/tmp/{sample}.amplicon.revcomp.fa',
        fq1='results/{experiment}/{sample}/tmp/{sample}_revcomp_R1_001.fastq.gz',
        fq2='results/{experiment}/{sample}/tmp/{sample}_revcomp_R2_001.fastq.gz'
    output:
        'results/{experiment}/{sample}/tmp/{sample}.amplicon.revcomp.fa.bwameth.c2t',
        'results/{experiment}/{sample}/tmp/{sample}.amplicon.revcomp.fa.bwameth.c2t.amb', 
        'results/{experiment}/{sample}/tmp/{sample}.amplicon.revcomp.fa.bwameth.c2t.sa' 
    conda:
        "envs/python3_v7.yaml"
    shell:
        'amplicon-smf/workflow/scripts/bwameth.py index {input.fa}'
        # 'ml python/2.7.13; set +u; source $HOME/bin/VENV/VIRTUALENV_2.7.13_SAMSTATS/bin/activate; set -u; $HOME/bin/copy_from_georgi/bwameth.py index {input.fa}'
        # '$HOME/bin/bwa-meth/bwameth.py index {input.fa}'

rule align_bwameth:
    input:
        amplicon='results/{experiment}/{sample}/tmp/{sample}.amplicon.revcomp.fa',
        amplicon_index='results/{experiment}/{sample}/tmp/{sample}.amplicon.revcomp.fa.bwameth.c2t',
        amplicon_index_tmp1='results/{experiment}/{sample}/tmp/{sample}.amplicon.revcomp.fa.bwameth.c2t.amb', 
        amplicon_index_tmp2='results/{experiment}/{sample}/tmp/{sample}.amplicon.revcomp.fa.bwameth.c2t.sa', 
        read1='results/{experiment}/{sample}/tmp/{sample}_revcomp_R1_001.fastq.gz',
        read2='results/{experiment}/{sample}/tmp/{sample}_revcomp_R2_001.fastq.gz'
        # read1='results/{sample}/tmp/{sample}.read1.revcomp.fastq.gz',
        # read2='results/{sample}/tmp/{sample}.read2.revcomp.fastq.gz'
    output:
        sam=temp('results/{experiment}/{sample}/tmp/{sample}.bwameth.sam')
    log:
        'results/{experiment}/{sample}/tmp/{sample}.bwameth.log'
    params:
        threads=config.get('threads', 1) # 1
    conda:
        "envs/python3_v7.yaml"
    shell:
        # '/oak/stanford/groups/wjg/bgrd/scripts_share/bwameth.py --threads {params.threads} --reference {input.amplicon} {input.read1} {input.read2} > {output.sam}'
        # '/oak/stanford/groups/wjg/bgrd/scripts_share/bwameth.py --threads {params.threads} --reference {input.amplicon} {input.read1} {input.read2} | grep -v "[W::sam_parse1]" > {output.sam}'
        'amplicon-smf/workflow/scripts/bwameth_nobenedit.py --threads {params.threads} --reference {input.amplicon} {input.read1} {input.read2} > {output.sam} 2> {log} || (echo "align_bwameth failed; stderr follows:" >&2 ; cat {log} >&2 ; exit 1)'
        # '$HOME/bin/copy_from_georgi/bwameth.py --threads {params.threads} --reference {input.amplicon} {input.read1} {input.read2} > {output.sam}'
        # 'ml python/2.7.13; set +u; source $HOME/bin/VENV/VIRTUALENV_2.7.13_SAMSTATS/bin/activate; set -u; $HOME/bin/copy_from_georgi/bwameth.py --threads {params.threads} --reference {input.amplicon} {input.read1} {input.read2} > {output.sam}'
        # '$HOME/bin/bwa-meth/bwameth.py --threads {params.threads} --reference {input.amplicon} {input.read1} {input.read2} > {output.sam}'

rule align_bwameth_all:
    input:
        amplicon='results/{experiment}/{sample}/tmp/{sample}.amplicon.revcomp.fa',
        amplicon_index='results/{experiment}/{sample}/tmp/{sample}.amplicon.revcomp.fa.bwameth.c2t',
        amplicon_index_tmp1='results/{experiment}/{sample}/tmp/{sample}.amplicon.revcomp.fa.bwameth.c2t.amb', 
        amplicon_index_tmp2='results/{experiment}/{sample}/tmp/{sample}.amplicon.revcomp.fa.bwameth.c2t.sa', 
        read1='results/{experiment}/{sample}/tmp/{sample}_revcomp_R1_001.fastq.gz',
        read2='results/{experiment}/{sample}/tmp/{sample}_revcomp_R2_001.fastq.gz'
        # read1='results/{sample}/tmp/{sample}.read1.revcomp.fastq.gz',
        # read2='results/{sample}/tmp/{sample}.read2.revcomp.fastq.gz'
    output:
        sam=temp('results/{experiment}/{sample}/tmp/{sample}.bwameth.all_alignments.sam')
    params:
        threads=config.get('threads', 1)
    conda:
        "envs/python3_v7.yaml"
    shell:
        #'ml python/2.7.13; set +u; source $HOME/bin/VENV/VIRTUALENV_2.7.13_SAMSTATS/bin/activate; set -u; $HOME/bin/copy_from_georgi/bwameth.py --threads {params.threads} --reference {input.amplicon} {input.read1} {input.read2} > {output.sam}'
        'amplicon-smf/workflow/scripts/bwameth.py --threads {params.threads} --reference {input.amplicon} {input.read1} {input.read2} > {output.sam}'
        # 'amplicon-smf/workflow/scripts/bwameth.py --threads {params.threads} --reference {input.amplicon} {input.read1} {input.read2} | grep -v "[W::sam_parse1]" > {output.sam}'
        # 'ml python/2.7.13; set +u; source $HOME/bin/VENV/VIRTUALENV_2.7.13_SAMSTATS/bin/activate; set -u; /oak/stanford/groups/wjg/bgrd/projects/smf/210817_smf_code_rewrite_reanalysis/210825_rerun_all_smf_snakemake/pipeline_test/bwameth.py --threads {params.threads} --reference {input.amplicon} {input.read1} {input.read2} > {output.sam}'
        # 'amplicon-smf/workflow/scripts/bwameth.py --threads {params.threads} --reference {input.amplicon} {input.read1} {input.read2} > {output.sam}'

# rule filter_sam_file:
#     input: 
#         'results/{experiment}/{sample}/tmp/{sample}.bwameth.all_alignments.sam'
#     output:
#         'results/{experiment}/{sample}/tmp/{sample}.bwameth.all_alignments.filtered.sam'
#     conda:
#         "envs/python3_v7.yaml"
#     shell:
#         'samtools view {input} | grep -v "[W::sam_parse1]" > {output}'

# need this to determine whether we should filter out the contigs like 0x-5x or do normal
def choose_bwameth_version(wildcards):
    if samplesheet.loc[wildcards.sample, 'filter_contigs']:
        # return 'results/{}/{}/tmp/{}.bwameth.all_alignments.sam'.format(wildcards.experiment, wildcards.sample, wildcards.sample)
        return 'results/{}/{}/tmp/{}.bwameth.all_alignments.sam'.format(wildcards.experiment, wildcards.sample, wildcards.sample)
    else:
        return 'results/{}/{}/tmp/{}.bwameth.sam'.format(wildcards.experiment, wildcards.sample, wildcards.sample)

rule correct_mismatched_amplicons:
    input:
        bam=choose_bwameth_version,
        fa='results/{experiment}/{sample}/tmp/{sample}.amplicon.revcomp.fa'
        # bam='results/{experiment}/{sample}/tmp/{sample}.bwameth.all_alignments.sam'
    output:
        bam=temp('results/{experiment}/{sample}/tmp/{sample}.bwameth.contig_filtered.sam'),
        out_stats='results/{experiment}/{sample}/stats/{sample}.bwameth.contig_filtered.stats.txt',
        failure_modes='results/{experiment}/{sample}/stats/{sample}.bwameth.contig_filtered.failure_modes.txt',
        problematic_reads=temp('results/{experiment}/{sample}/tmp/{sample}.bwameth.problematic_reads.sam')
    params:
        read1_thresh=lambda wildcards: int(samplesheet.loc[wildcards.sample,'read1_length'] * config['alignment_length_fraction']), #  config['read1_length_override'].get(wildcards.sample, config['read1_length_threshold']) if 'read1_length_override' in config.keys() else config['read1_length_threshold'],
        read2_thresh=lambda wildcards: int(samplesheet.loc[wildcards.sample,'read2_length'] * config['alignment_length_fraction']), # config['read2_length_override'].get(wildcards.sample, config['read2_length_threshold']) if 'read2_length_override' in config.keys() else config['read2_length_threshold'],
        as_thresh=config['alignment_score_fraction'],
        ignore_bounds=lambda wildcards: '--ignore_bounds' if samplesheet.loc[wildcards.sample,'ignore_bounds'] else '', # wildcards.sample in config['ignore_bounds'] else ''
        write_problematic_reads=lambda wildcards: '--write_problematic_reads' if True else '' 
    conda:
        "envs/python3_v7.yaml"
    shell:
        'python amplicon-smf/workflow/scripts/filter_bam_by_matching_contigs2.py --bam {input.bam} --out {output.bam} --out_stats {output.out_stats} --failure_modes {output.failure_modes} --problem {output.problematic_reads} --amplicon {input.fa} --min_len_threshold1 {params.read1_thresh} --min_len_threshold2 {params.read2_thresh} --min_as_frac {params.as_thresh} {params.ignore_bounds} {params.write_problematic_reads}'

rule plot_mapped_vs_unmapped_reads:
    input:
        # ot_bam='results/{experiment}/{sample}/tmp/{sample}.bwameth.contig_filtered.sam',
        # problematic_reads='results/{experiment}/{sample}/tmp/{sample}.bwameth.problematic_reads.sam'
        'results/{experiment}/{sample}/stats/{sample}.bwameth.contig_filtered.stats.txt'
    output:
        'results/{experiment}/plots/{sample}.wasted_reads.pdf'
    conda:
        "envs/python3_v7.yaml"
    shell:
        # 'python amplicon-smf/workflow/scripts/count_excluded_reads.py --retained {input.ot_bam} --removed {input.problematic_reads} --plot {output}'
        'python amplicon-smf/workflow/scripts/count_excluded_reads.py --stats {input} --plot {output}'

# # need this to determine whether we should filter out the contigs like 0x-5x or do normal
# def choose_contig_filtered_bam(wildcards):
#     if samplesheet.loc[wildcards.sample, 'filter_contigs']:
#         return 'results/{}/{}/tmp/{}.bwameth.contig_filtered.sam'.format(wildcards.experiment, wildcards.sample, wildcards.sample)
#     else:
#         return 'results/{}/{}/tmp/{}.bwameth.sam'.format(wildcards.experiment, wildcards.sample, wildcards.sample)

rule sam_to_bam:
    input:
        # sam='results/{experiment}/{sample}/tmp/{sample}.bwameth.sam',
        # sam='results/{experiment}/{sample}/tmp/{sample}.bwameth.contigs.sam',
        sam='results/{experiment}/{sample}/tmp/{sample}.bwameth.contig_filtered.sam',
        # sam=choose_contig_filtered_bam,
        fa='results/{experiment}/{sample}/tmp/{sample}.amplicon.revcomp.fa'
    output:
        'results/{experiment}/{sample}/{sample}.bwameth.bam'
    conda:
        "envs/python3_v7.yaml"
    shell:
        'samtools view -F 1804 -q 30 -bT {input.fa} {input.sam} > {output}'

rule sort_index_bam:
    input:
        'results/{experiment}/{sample}/{sample}.{extra}.bam'
    output:
        'results/{experiment}/{sample}/{sample}.{extra}.sorted.bam'
    conda:
        "envs/python3_v7.yaml"
    shell:
        '''
        samtools sort {input} > {output}
        samtools index {output}
        '''

rule filter_uncoverted:
    input:
        bam='results/{experiment}/{sample}/{sample}.bwameth.sorted.bam',
        fa='results/{experiment}/{sample}/tmp/{sample}.amplicon.revcomp.fa'
    output:
        bam='results/{experiment}/{sample}/{sample}.bwameth.filtered.bam',
        plot='results/{experiment}/plots/{sample}.nonconverted_reads.png'
    params:
        cfrac=lambda wildcards: 0.001 if samplesheet.loc[wildcards.sample, 'deaminase'] else config.get('unconverted_frac', 0.85) 
    conda:
        "envs/python3_v7.yaml"
    shell:
        'python amplicon-smf/workflow/scripts/mark-nonconverted-reads-and-plot.py --reference {input.fa} --bam {input.bam} --out {output.bam} --plot {output.plot} --c_frac {params.cfrac}'

# rule methylation_percentage:
#     input:
#     output:
#     shell:

rule run_methyldackel:
    input:
        bam='results/{experiment}/{sample}/{sample}.bwameth.filtered.sorted.bam',
        fa='results/{experiment}/{sample}/tmp/{sample}.amplicon.revcomp.fa'
    output:
        'results/{experiment}/{sample}/{sample}.bwameth.filtered.sorted_CpG.bedGraph'
    conda:
        "envs/python3_v7.yaml"
    shell:
        # 'MethylDackel extract --CHG --CHH --keepSingleton --keepDiscordant {input.fa} {input.bam}'
        # 'MethylDackel extract --CHG --CHH --ignoreFlags 3584 {input.fa} {input.bam}'
        # 'MethylDackel extract --CHG --CHH  --keepSingleton --ignoreFlags 3584 {input.fa} {input.bam}'
        'MethylDackel extract --CHG --CHH {input.fa} {input.bam}'
        # figure this out


rule plot_bulk_methylation:
    input:
        bg='results/{experiment}/{sample}/{sample}.bwameth.filtered.sorted_CpG.bedGraph',
        fa='results/{experiment}/{sample}/tmp/{sample}.amplicon.revcomp.fa'
    output:
        'results/{experiment}/plots/{sample}.bulk_plots.pdf'
    params:
        prefix='results/{experiment}/{sample}/{sample}.bwameth.filtered.sorted',
        cpg=lambda wildcards: '--include_cpg' if samplesheet.loc[wildcards.sample, 'include_cpg'] else '',
        no_endog_meth=lambda wildcards: '--no_endog_meth' if samplesheet.loc[wildcards.sample, 'no_endog_meth'] else '',
        deaminase=lambda wildcards: '--deaminase' if samplesheet.loc[wildcards.sample, 'deaminase'] else ''
    conda:
        "envs/python3_v7.yaml"
    shell:
        'python amplicon-smf/workflow/scripts/plot_bulk_methylation_signal.py --input {params.prefix} --amplicon {input.fa} --plot {output} {params.cpg} {params.no_endog_meth} {params.deaminase}'

rule amplicon_fa_to_peak_bed:
    input:
        'results/{experiment}/{sample}/tmp/{sample}.amplicon.revcomp.fa'
    output:
        'results/{experiment}/{sample}/tmp/{sample}.amplicon.revcomp.peaklist.bed'
    conda:
        "envs/python3_v7.yaml"
    shell:
        'python amplicon-smf/workflow/scripts/convert_amplicon_fa_to_peaklist.py {input} {output}'

def get_c_type(wildcards):
    if samplesheet.loc[wildcards.sample,'deaminase']:
        return 'allC'
    else:
        if samplesheet.loc[wildcards.sample,'include_cpg']:
            return 'both_dimers'
        else:
            return 'GC'

rule join_reads_and_first_cluster:
    input:
        bam='results/{experiment}/{sample}/{sample}.bwameth.filtered.sorted.bam',
        fa='results/{experiment}/{sample}/tmp/{sample}.amplicon.revcomp.fa',
        peaks='results/{experiment}/{sample}/tmp/{sample}.amplicon.revcomp.peaklist.bed'
    output:
        'results/{experiment}/{sample}/{sample}.amplicon_stats.txt'
    params:
        prefix='results/{experiment}/{sample}/matrices/{sample}',
        matdir='results/{experiment}/{sample}/matrices',
        subset=1000,
        # ctype=lambda wildcards: 'both' if samplesheet.loc[wildcards.sample,'include_cpg'] else 'GC',
        ctype=get_c_type,
        no_endog_meth=lambda wildcards: '-noEndogenousMethylation' if samplesheet.loc[wildcards.sample, 'no_endog_meth'] else ''
    conda:
        "envs/python3_v7.yaml"
    shell:
        'mkdir -p {params.matdir}; python amplicon-smf/workflow/scripts/dSMF_footprints_clustering_py3.py {input.bam} {input.fa} {params.ctype} {input.peaks} 0 1 2 3 {params.prefix} {output} -label 0 -unstranded -subset {params.subset} {params.no_endog_meth} -cluster -heatmap'

rule plot_bulk_methylation2:
    input:
        fa='results/{experiment}/{sample}/tmp/{sample}.amplicon.revcomp.fa',
        stats='results/{experiment}/{sample}/{sample}.amplicon_stats.txt'
    output:
        'results/{experiment}/plots/{sample}.bulk_plots_from_matrices.pdf'
    params:
        matrix_path='results/{experiment}/{sample}/matrices/{sample}'
    conda:
        "envs/python3_v7.yaml"
    shell:
        'python amplicon-smf/workflow/scripts/plot_bulk_methylation_from_matrices.py --input {params.matrix_path} --amplicon {input.fa} --plot {output}'

rule plot_nuc_qc:
    input:
        fa='results/{experiment}/{sample}/tmp/{sample}.amplicon.revcomp.fa',
        stats='results/{experiment}/{sample}/{sample}.amplicon_stats.txt'
    output:
        plot='results/{experiment}/plots/{sample}.nuc_len_qc_plots.pdf',
        stats='results/{experiment}/{sample}/stats/{sample}.nuc_len_qc.stats.txt'
    params:
        matrix_path='results/{experiment}/{sample}/matrices/{sample}'
    conda:
        "envs/python3_v7.yaml"
    shell:
        'python amplicon-smf/workflow/scripts/plot_nucleosome_qc.py --input {params.matrix_path} --amplicon {input.fa} --plot {output.plot} --results {output.stats} --gmm'

# rule collate_stats:

# # emul specific analyses (glue vertically, as in combine stats from a bunch of separate places into one)
# rule combine_emul_stats:
#     input:
#         dup_stats_un="{sample}/stats/{sample}_uncollapsed.duplication_stats.txt",
#         dup_stats="{sample}/stats/{sample}_collapsed.duplication_stats.txt",
#         dup_stats_filt="{sample}/stats/{sample}_collapsed_filtered.duplication_stats.txt",
#         count_stats_un="{sample}/stats/{sample}_uncollapsed.out_stats.txt",
#         count_stats="{sample}/stats/{sample}_collapsed.out_stats.txt",
#         count_stats_filt="{sample}/stats/{sample}_collapsed_filtered.out_stats.txt",
#         collapse_stats="{sample}/stats/{sample}.collapse_barcodes.stats.txt",
#         collapsed_reads="{sample}/stats/{sample}_emul_extracted_collapsed.stats.txt",
#         uncollapsed_reads="{sample}/stats/{sample}_emul_extracted_uncollapsed.stats.txt",
#         collapsed_filt_reads="{sample}/stats/{sample}_emul_extracted_collapsed_filtered.stats.txt",
#         filtered_reads="{sample}/stats/{sample}_filtered.stats.txt",
#         input_reads="{sample}/stats/{sample}_downsampled.stats.txt",
#         starting_cells="{sample}/stats/{sample}_starting_cells.stats.txt"
#     output:
#         "{sample}/{sample}.emul_stats.txt"
#     shell:
#         "python {mip_directory}/combine_stats.py -i {input.dup_stats} {input.dup_stats_un} {input.dup_stats_filt} {input.count_stats} {input.count_stats_un} {input.count_stats_filt} \
#             {input.collapse_stats} {input.collapsed_reads} {input.uncollapsed_reads} {input.collapsed_filt_reads} {input.filtered_reads} {input.input_reads} \
#             {input.starting_cells} -o {output} -e .emul"

# rule emul_summary:
#     input:
#         counts=["{name}-{id}/{name}-{id}_collapsed_emul.sum_genes.txt".format(name=row.Name, id=row.ID) for row in sample_sheet.loc[sample_sheet.Read2 != 'Bulk'].itertuples()],
#         stats=["{name}-{id}/{name}-{id}.emul_stats.txt".format(name=row.Name, id=row.ID) for row in sample_sheet.loc[sample_sheet.Read2 != 'Bulk'].itertuples()]
#     output:
#         counts_summary='emul_count_summary.txt',
#         stats_summary='emul_stats_summary.txt'
#     params:
#         emul_summarize=path.join(mip_directory, 'emul_summarize.py')
#     shell:
#         "python {params.emul_summarize} --output_stats {output.stats_summary} --input_stats {input.stats} \
#             --output_counts {output.counts_summary} --input_counts {input.counts}"


