configfile: "config.yaml"

rule all: 
    input:
        expand("{sample}.narrowPeak", sample=config['SAMPLES'])
 
rule trim: 
    input: 
       r1 = "{sample}.r_1.fq.gz",
       r2 = "{sample}.r_2.fq.gz"
    output: 
      val1 = "galore/{sample}.r_1_val_1.fq.gz",
      val2 = "galore/{sample}.r_2_val_2.fq.gz"
    shell: 
        """
         trim_galore --gzip --retain_unpaired --trim1 --fastqc --fastqc_args "--outdir fastqc" -o galore --paired {input.r1} {input.r2}
        """ 

rule tosam:
    input:
        genome = config['GENOME'],
        r1 = "galore/{sample}.r_1_val_1.fq.gz",
        r2 = "galore/{sample}.r_2_val_2.fq.gz"
    output:
        '{sample}.sam'
    shell:
        "bowtie2 -x {input.genome} -1 {input.r1} -2 {input.r2} -S {output}"

rule tobam:
      input:
          "{sample}.sam"
      output:
          "{sample}.bowtie2.bam"
      shell:
          "samtools view {input[0]} -S -b > {output[0]}"

rule sort:
    input: 
       "{sample}.bowtie2.bam"
    output:
       "{sample}.sorted.bam" 
    params: 
        "{sample}.tmp.sorted"
    shell: 
       "samtools sort -T {params}  -o {output} {input};"

rule peak_call: 
    input:
       "{sample}.sorted.bam" 
    params: 
      "{sample}.log",
      "{sample}.err",
      expand("{chr}", chr=config['CHR'])
    output: 
       "{sample}.narrowPeak"
    shell: 
        "Genrich -t {input} -o {output}  -r -j -v -e {params[2]} > {params[0]} 2> {params[1]}"  
