rule fastqc:
    input:
        r1 = lambda wildcards: FILES[wildcards.sample]['R1'],
        r2 = lambda wildcards: FILES[wildcards.sample]['R2'],
    output:
        r1 = temp(join(OUT_DIR, 'qc', 'fastqc', '{sample}.R1_fastqc.html')),
        r2 = temp(join(OUT_DIR, 'qc', 'fastqc', '{sample}.R2_fastqc.html')),
    params:
        out_dir = join(OUT_DIR, 'qc', 'fastqc')
    log:
        join(OUT_DIR, 'logs', 'fastqc', '{sample}.fastQC.log')
    benchmark:
        join(OUT_DIR, 'benchmarks', 'fastqc', '{sample}.fastQC.benchmark.tsv')
    message:
        """--- Checking read quality of sample {wildcards.sample} with FastQC """
    conda:
        "../envs/rnaseq_yablab.yml"
    shell:
        """
        fastqc -o {params.out_dir} -t 4 -q {input.r1} {input.r2} > {log} 2>&1
        """

rule multiQC:
    input:
        r1 = expand(join(OUT_DIR, 'qc', 'fastqc', '{sample}.R1_fastqc.html'), sample = SAMPLES),
        r2 = expand(join(OUT_DIR, 'qc', 'fastqc', '{sample}.R2_fastqc.html'), sample = SAMPLES),
        qc = expand(join(OUT_DIR, 'trimmed_reads', '{sample}.qc.txt'), sample = SAMPLES),
        featureCounts = join(OUT_DIR, 'logs', 'abundances', 'featureCounts_fC.log'),
    output:
        join(OUT_DIR, 'qc', 'multiqc_report.html')
    log:
        join(OUT_DIR, 'qc', 'multiQC.log')
    benchmark:
        join(OUT_DIR, 'qc', 'multiQC.benchmark.tsv')
    params:
        summary_file = join(OUT_DIR, 'qc', 'summary_files.txt'),
        fastqc_dir = join(OUT_DIR, 'qc', 'fastqc', '*.fastqc.zip'),
        trim_dir = join(OUT_DIR, 'trimmed_reads', '*.qc.txt'),
        out_dir = join(OUT_DIR, 'qc')
    conda:
        "../envs/rnaseq_yablab.yml"
    message:
        """--- Running MultiQC """
    shell:
        """
        ls {params.fastqc_dir} > {params.summary_file}
        ls {params.trim_dir} >> {params.summary_file}
        echo "{input.featureCounts}" >> {params.summary_file}
        multiqc -f -o {params.out_dir} -d -dd 3 -l {params.summary_file} > {log} 2>&1
        """