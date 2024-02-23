

# run snakemake
snakemake -s amplicon-smf/workflow/Snakefile -c1 -k -w 15 --configfile amplicon-smf/example_files/example_config.yaml --use-conda


