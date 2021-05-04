Author: Sherine Awad 

A snakefile pipeline for ATAC-Seq.

Change the config.yaml file appropriately according to your data. 
For example, change workdir, reference genome, samples prefixes in replicates, etc. 

Then run: snakemake -jnumber_of_cores, for example for 5 cores use:

    snakemake -j5 

and for a dry run use: 

    snakemake -j1 -n 


and to print the commands in a dry run use:

    snakemake -j1 -n -p 

