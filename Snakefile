configfile: "config.yaml"
workdir: config['WORKDIR']

import glob
from pathlib import Path
REPLICATE = [Path(Path(Path(x).stem).stem).stem for x in glob.glob("*r_1.fq.gz")]

rule all: 
    input:
        expand("{sample}.sam", sample = REPLICATE), 
        expand("{sample}.bam", sample = REPLICATE), 
        expand("{sample}.sorted.bam", sample = REPLICATE),    
        expand("{replicate}.narrowPeak", replicate = config['REPLICATE_NAME']) 

rule trim: 
    input: 
       r1 = "{sample}.r_1.fq.gz",
       r2 = "{sample}.r_2.fq.gz"
    output: 
      val1 = "galore/{sample}.r_1_val_1.fq.gz",
      val2 = "galore/{sample}.r_2_val_2.fq.gz"
    log: 
        "{sample}.trim.log"
    shell: 
         """
         trim_galore --gzip --retain_unpaired --trim1 --fastqc --fastqc_args "--outdir fastqc" -o galore --paired {input.r1} {input.r2}
         """ 

rule tosam:
    input:
        r1 = "galore/{sample}.r_1_val_1.fq.gz",
        r2 = "galore/{sample}.r_2_val_2.fq.gz"
    params:
       genome = config['GENOME']
    output:
        "{sample}.sam"
    log: 
        "{sample}.sam.log"
    shell:
        "bowtie2 -x {params} -1 {input.r1} -2 {input.r2} -S {output}"

rule tobam:
      input:
          "{sample}.sam"
      output:
          "{sample}.bam"
      log: 
          "{sample}.bam.log"
      shell:
          "samtools view {input[0]} -S -b > {output[0]}"

rule sort:
    input: 
       "{sample}.bam"
    output:
       "{sample}.sorted.bam" 
    params: 
        "{sample}.tmp.sorted"
    log: 
        "{sample}.sorted.log" 
    shell: 
       "samtools sort -T {params} -n -o {output} {input}"

rule peak_call:
    input: 
        expand("{sample}.sorted.bam", sample = REPLICATE) 
    params:
       lambda w: ",".join(expand("{sample}.sorted.bam", sample = REPLICATE)), 
       expand("{replicate}.log", replicate = config['REPLICATE_NAME']),
       expand("{replicate}.err", replicate = config['REPLICATE_NAME']),
       expand("{chr}", chr=config['CHR'])
    output: 
       expand("{replicate}.narrowPeak", replicate = config['REPLICATE_NAME']) 
    shell:
       "Genrich -t {params[0]} -o {output}  -r -j -v -e {params[3]} > {params[1]} 2> {params[2]}"  
