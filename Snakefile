configfile: "config.yaml"
workdir: config['WORKDIR']

rule all: 
     input:
        expand("galore/{sample}.r_1_val_1.fq.gz", sample = config['REPLICATE']), 
        expand("{sample}.sam", sample = config['REPLICATE']), 
        expand("{sample}.bam", sample = config['REPLICATE']), 
        expand("{sample}.sorted.bam", sample = config['REPLICATE']),    
        expand("{replicate}.narrowPeak", replicate = config['REPLICATE_NAME']) 

rule trim: 
    input: 
       r1 = "{sample}.r_1.fq.gz",
       r2 = "{sample}.r_2.fq.gz"
    output: 
      val1 = "galore/{sample}.r_1_val_1.fq.gz",
      val2 = "galore/{sample}.r_2_val_2.fq.gz"
    shell: 
         "trim_galore --gzip --retain_unpaired --trim1 --fastqc --fastqc_args "--outdir fastqc" -o galore --paired {input.r1} {input.r2}" 

rule tosam:
    input:
        r1 = "galore/{sample}.r_1_val_1.fq.gz",
        r2 = "galore/{sample}.r_2_val_2.fq.gz"
    params:
       genome = config['GENOME']
    output:
        '{sample}.sam'
    shell:
        "bowtie2 -x {params} -1 {input.r1} -2 {input.r2} -S {output}"

rule tobam:
      input:
          "{sample}.sam"
      output:
          "{sample}.bam"
      shell:
          "samtools view {input[0]} -S -b > {output[0]}"

rule sort:
    input: 
       "{sample}.bam"
    output:
       "{sample}.sorted.bam" 
    params: 
        "{sample}.tmp.sorted"
    shell: 
       "samtools sort -T {params} -n -o {output} {input}"

rule peak_call: 
    params:
       lambda w: ",".join(expand("{sample}.sorted.bam", sample = config['REPLICATE'])), 
       expand("{replicate}.log", replicate = config['REPLICATE_NAME']),
       expand("{replicate}.err", replicate = config['REPLICATE_NAME']),
       expand("{chr}", chr=config['CHR'])
    output: 
       expand("{replicate}.narrowPeak", replicate = config['REPLICATE_NAME']) 
    shell:
       "Genrich -t {params[0]} -o {output}  -r -j -v -e {params[3]} > {params[1]} 2> {params[2]}"  
