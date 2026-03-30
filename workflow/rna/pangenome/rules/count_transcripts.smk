rule build_index:
    input: 
        graph_gbz = config["reference_gbz"],
        gtf = GTF
    output: 
        xg = temp(join(OUT_DIR, "pantranscriptome", "dame_graph.xg")),
        gcsa = temp(join(OUT_DIR, "pantranscriptome", "dame_graph.gcsa"))
    params:
        out_dir = join(OUT_DIR, "pantranscriptome", "dame_graph"),
        tmp_dir = join(OUT_DIR, "temp"),
    log:
        join(OUT_DIR, "logs", "index.txt")
    resources:
        mem="1000GB"
    threads:
        32
    conda:
        "../envs/cactus-minigraph.yml"
    message:
        """--- Building RNA index """
    shell:
        """
        vg autoindex -w mpmap --threads {threads} --target-mem 1000 -p {params.out_dir} -T {params.tmp_dir} -f exon --gbz {input.graph_gbz} -x {input.gtf} 2> {log}
        """

rule cutadapt:
    input:
        lambda wildcards: FILES[wildcards.sample]['R1'],
        lambda wildcards: FILES[wildcards.sample]['R2']
    output:
        fastq1 = temp(join(OUT_DIR, 'trimmed_reads', '{sample}.1.fq.gz')),
        fastq2 = temp(join(OUT_DIR, 'trimmed_reads', '{sample}.2.fq.gz')),
        qc = temp(join(OUT_DIR, 'trimmed_reads', '{sample}.qc.txt'))
    params:
        # https://cutadapt.readthedocs.io/en/stable/guide.html#adapter-types
        adapters="-a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -a CTGTCTCTTATACACATCT -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -A CTGTCTCTTATACACATCT",
        # https://cutadapt.readthedocs.io/en/stable/guide.html#
        extra="--minimum-length 30 -q 20"
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
        "v4.4.0/bio/cutadapt/pe"

rule mpmap_bam:
    input:
        fastq1 = join(OUT_DIR, 'trimmed_reads', '{sample}.1.fq.gz'),
        fastq2 = join(OUT_DIR, 'trimmed_reads', '{sample}.2.fq.gz'),
        xg = join(OUT_DIR, "pantranscriptome", "dame_graph.xg"),
        gcsa = join(OUT_DIR, "pantranscriptome", "dame_graph.gcsa")
    output:
        mapped_bam = temp(join(OUT_DIR, "abundances", "{sample}"))
    params:
        paths = temp(join(OUT_DIR, "abundances", "dame_paths.txt"))
    message:
        """--- Align {wildcards.sample} with VG """
    conda:
        "../envs/cactus-minigraph.yml"
    shell:
        """
        mkdir -p {output.mapped_bam}
        vg paths -L -x {input.xg} > {params.paths}
        vg mpmap -x {input.xg} -g {input.gcsa} -n RNA -f {input.fastq1} -f {input.fastq2} -F BAM --ref-paths {params.paths} > {output.mapped_bam}
        """

rule featureCounts:
    input:
        bams = expand(join(OUT_DIR, "abundances", "{sample}"), sample=SAMPLES)
    output:
        countsSt = join(OUT_DIR, 'abundances', 'gene_counts.txt'),
        countsStMatrix = join(OUT_DIR, 'abundances', 'gene_counts_matrix.txt')
    params:
        gtf = GTF
    log:
        join(OUT_DIR, 'logs', 'abundances', 'featureCounts_fC.log')
    benchmark:
        join(OUT_DIR, 'benchmarks', 'abundances', 'featureCounts_fC.benchmark.tsv')
    threads:
        8
    resources:
        mem_mb=16000
    conda:
        "../envs/rnaseq_yablab.yml"
    message:
        """--- Outputting gene counts with featureCounts for genome mapping. """
    shell:
        """
        featureCounts -a {params.gtf} -p -T {threads} -g "gene_name" -F GTF -o {output.countsSt} {input.bams} > {log} 2>&1
        cat {output.countsSt} | cut -f 1,7- | sed 1d > {output.countsStMatrix}
        """