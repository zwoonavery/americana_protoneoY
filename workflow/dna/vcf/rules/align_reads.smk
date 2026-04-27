rule cutadapt_trim:
    input:
        lambda wildcards: FILES[wildcards.sample]['R1'],
        lambda wildcards: FILES[wildcards.sample]['R2']
    output:
        fastq1 = temp(join(OUT_DIR, 'trimmed_reads', '{sample}.R1.fastq.gz')),
        fastq2 = temp(join(OUT_DIR, 'trimmed_reads', '{sample}.R2.fastq.gz')),
        qc = join(OUT_DIR, 'qc', '{sample}', '{sample}.qc.txt')
    params:
        # https://cutadapt.readthedocs.io/en/stable/guide.html#adapter-types
        adapters="-a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -a CTGTCTCTTATACACATCT -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -A CTGTCTCTTATACACATCT",
        # https://cutadapt.readthedocs.io/en/stable/guide.html#
        extra="--minimum-length 50 -q 20"
    log:
        join(OUT_DIR, 'logs', 'cutadapt', '{sample}.cutadapt.log')
    benchmark:
        join(OUT_DIR, 'benchmarks', 'cutadapt', '{sample}.cutadapt.benchmark.tsv')
    message:
        """--- Trimming adaptors for sample {wildcards.sample}."""
    threads:
        4
    resources:
        mem_mb=8000
    wrapper:
        "v3.10.2/bio/cutadapt/pe"

rule vg_pack:
    input:
        fastq1 = join(OUT_DIR, 'trimmed_reads', '{sample}.R1.fastq.gz'),
        fastq2 = join(OUT_DIR, 'trimmed_reads', '{sample}.R2.fastq.gz')
    output:
        vg_gam = temp(join(OUT_DIR, 'trimmed_reads', '{sample}.gam')),
    log: 
        join(OUT_DIR, 'logs', 'align', '{sample}_pack.log')
    params:
        giraffe_gbz = expand(join(GIRAFFE_FOLDER, '{prefix}.giraffe.gbz'), prefix = PREFIX),
        giraffe_dist = expand(join(GIRAFFE_FOLDER, '{prefix}.dist'), prefix = PREFIX),
        giraffe_zip = expand(join(GIRAFFE_FOLDER, '{prefix}.shortread.zipcodes'), prefix = PREFIX),
        giraffe_min = expand(join(GIRAFFE_FOLDER, '{prefix}.shortread.withzip.min'), prefix = PREFIX)
    threads:
        32
    message:
        """--- Package {wildcards.sample} reads into GAM """
    conda:
        '../envs/vg.yml'
    shell:
        """
        vg giraffe -p -Z {params.giraffe_gbz} -d {params.giraffe_dist} -z {params.giraffe_zip} -m {params.giraffe_min} \
        -f {input.fastq1} -f {input.fastq2} > {output.vg_gam} 2> {log}
        """

rule vg_align:
    input:
        vg_gam = join(OUT_DIR, 'trimmed_reads', '{sample}.gam'),
        giraffe_gbz = expand(join(GIRAFFE_FOLDER, '{prefix}.giraffe.gbz'), prefix = PREFIX)
    output:
        vg_bam = temp(join(OUT_DIR, 'alignment', '{sample}_alignment.bam'))
    log: 
        join(OUT_DIR, 'logs', 'align', '{sample}_align.log')
    threads:
        32
    message:
        """--- Align {wildcards.sample} to pangenome """
    conda:
        '../envs/vg.yml'
    shell:
        """
        vg giraffe -Z {input.giraffe_gbz} -G {input.vg_gam} -t {threads} -o BAM > {output.vg_bam} 2> {log}
        """

rule sort_bam:
    input:
        vg_bam = join(OUT_DIR, 'alignment', '{sample}_alignment.bam')
    output:
        vg_sort = temp(join(OUT_DIR, 'alignment', '{sample}_sorted.bam')),
    threads:
        8
    resources:
        mem_mb=12000
    log:
        join(OUT_DIR, 'logs', 'align', '{sample}.sort.log')
    benchmark:
        join(OUT_DIR, 'benchmarks', 'align', '{sample}.sort.benchmark.tsv')
    message:
        """--- Sorting alignment for "{wildcards.sample}" with SamTools."""
    conda:
        '../envs/vcf.yml'
    shell:
        """
        samtools sort -@ {threads} -o {output.vg_sort} {input.vg_bam} 
        samtools index -@ {threads} {output.vg_sort}
        """

rule bamMarkDuplicates:
    input:
        vg_bam = join(OUT_DIR, 'alignment', '{sample}_alignment.bam')
    output:
        vg_mark = temp(join(OUT_DIR, 'alignment', '{sample}_dupMark.bam')),
    threads:
        8
    resources:
        mem_mb=12000
    log:
        join(OUT_DIR, 'logs', 'align', '{sample}.dupmark.log')
    benchmark:
        join(OUT_DIR, 'benchmarks', 'align', '{sample}.dupmark.benchmark.tsv')
    params:
        vg_sort = temp(join(OUT_DIR, 'alignment', '{sample}_sorted.bam'))
    message:
        """--- Marking duplicates for "{wildcards.sample}" with Picard."""
    conda:
        '../envs/vcf.yml'
    shell:
        """
        samtools sort -@ {threads} -o {params.vg_sort} {input.vg_bam} 
        picard MarkDuplicates -I {params.vg_sort} -O {output.vg_mark} -M {log} -VALIDATION_STRINGENCY LENIENT -ASSUME_SORTED true -REMOVE_DUPLICATES false
        samtools index -@ {threads} {output.vg_mark}
        """
