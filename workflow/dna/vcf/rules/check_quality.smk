rule check_reads:
    input:
        r1 = lambda wildcards: FILES[wildcards.sample]['R1'],
        r2 = lambda wildcards: FILES[wildcards.sample]['R2'],
    output:
        r1 = join(OUT_DIR, 'qc', 'fastqc' '{sample}.R1_fastqc.html'),
        r2 = join(OUT_DIR, 'qc', 'fastqc', '{sample}.R2_fastqc.html'),
    params:
        fastqc_folder = join(OUT_DIR, 'fastqc')
    log:
        join(OUT_DIR, 'logs', 'fastqc', '{sample}.fastQC.log')
    benchmark:
        join(OUT_DIR, 'benchmarks', 'fastqc', '{sample}.fastQC.benchmark.tsv')
    message:
        """--- Checking read quality of sample "{wildcards.sample}" with FastQC """
    conda:
        '../envs/vcf.yml'
    shell:
        """
        fastqc -o {params.fastqc_folder} -t 4 -q {input.r1} {input.r2} > {log} 2>&1
        """

rule check_alignment:
    input:
        vg_bam = join(OUT_DIR, 'alignment', '{sample}_alignment.bam')
    output:
        bamqc = join(OUT_DIR, 'qc', '{sample}', 'bamqc', 'qualimapReport.html')
    params:
        sort_bam = temp(join(OUT_DIR, 'qc', '{sample}_alignment_sorted.bam')),
        bamqc_folder = join(OUT_DIR, 'qc', '{sample}', 'bamqc')
    log:
        bamqc = join(OUT_DIR, 'logs', 'align', '{sample}.bamqc.log')
    benchmark:
        join(OUT_DIR, 'benchmarks', 'align', '{sample}.benchmark.tsv')
    message:
        """--- Running BAMQC for sample "{wildcards.sample}"."""
    threads:
        8
    resources:
        mem_mb=32000
    conda:
        '../envs/vcf.yml'
    shell:
        """
        samtools sort -@ 8 -o {params.sort_bam} {input.vg_bam}
        qualimap bamqc -bam {params.sort_bam} -c -outdir {params.bamqc_folder} --java-mem-size=32G -nt {threads} 2> {log.bamqc}
        """

rule vg_novelty:
    input: 
        vg_gam = join(OUT_DIR, 'trimmed_reads', '{sample}.gam'),
        full_xg = expand(join(OUT_DIR, 'calling', 'vg', '{prefix}.xg'), prefix = PREFIX),
        snarls = expand(join(OUT_DIR, 'calling', 'vg', '{prefix}.snarls'), prefix = PREFIX),
        gff = config['REFERENCE_GTF'],
    output:
        annotated_novelty = join(OUT_DIR, 'calling', 'vg', '{sample}.annotated.tsv')
    params:
        annotated_novelty = temp(join(OUT_DIR, 'calling', 'vg', '{sample}.annotated_raw.tsv'))
    threads:
        16
    message:
        """--- Identify novel annotations in '{wildcards.sample}' with VG. """
    conda:
        '../envs/vg.yml'
    shell:
        """
        vg annotate -x {input.full_xg} -f {input.gff} -s {input.snarls} -a {input.vg_gam} -n -t {threads} > {params.annotated_novelty}
        sed '1d' {params.annotated_novelty} | awk '{{print $0 "\t{wildcards.sample}"}}' > {output.annotated_novelty}
        """

rule vg_novelty_concatenate:
    input:
        annotated_novelty = expand(join(OUT_DIR, 'calling', 'vg', '{sample}.annotated.tsv'), sample = SAMPLES)
    output:
        novelty_cat = join(OUT_DIR, 'calling', 'vg_novelty.tsv')
    message:
        """--- Concatenate sample novelty results. """
    shell:
        """
        cat {input.annotated_novelty} > {output.novelty_cat}
        sed -i '1i name\tlength.bp\tunaligned.bp\tknown.nodes\tknown.bp\tnovel.nodes\tnovel.bp\tsample' {output.novelty_cat}
        """

rule check_gatk:
    input:
        gatk_vcf = join(OUT_DIR, 'calling', 'gatk', '{sample}_gatk.vcf'),
        dna = config['REFERENCE_FASTA']
    output:
        evalGrp = join(OUT_DIR, 'qc', 'gatk', '{sample}.evalGrp')
    log:
        join(OUT_DIR, 'logs', 'gatk', '{sample}.VariantEval.log')
    benchmark:
        join(OUT_DIR, 'benchmarks', 'gatk', '{sample}.VariantEval.benchmark.tsv')
    message:
        """--- Evaluating variant calls with GATK for sample {wildcards.sample}."""
    threads:
        4
    resources:
        mem_mb=8000
    conda:
        '../envs/vcf.yml'
    shell:
        """
        gatk VariantEval -R {input.dna} -O {output.evalGrp} --eval {input.gatk_vcf}
        """

rule multiQC:
    input:
        fastqc1 = expand(join(OUT_DIR, 'qc', 'fastqc' '{sample}.R1_fastqc.html'), sample = SAMPLES),
        fastqc2 = expand(join(OUT_DIR, 'qc', 'fastqc', '{sample}.R2_fastqc.html'), sample = SAMPLES),
        qc = expand(join(OUT_DIR, 'qc', '{sample}', '{sample}.qc.txt'), sample = SAMPLES),
        bamqc = expand(join(OUT_DIR, 'qc', '{sample}', 'bamqc', 'qualimapReport.html'), sample = SAMPLES),
        # evalGrp = expand(join(OUT_DIR, 'qc', 'gatk', '{sample}.evalGrp'), sample = SAMPLES)
    output:
        multiqc = join(OUT_DIR, 'multiqc', 'multiqc_report.html')
    params:
        fastqc_zip = expand(join(OUT_DIR, 'qc', 'fastqc' '*fastqc.zip'), sample = SAMPLES),
        multiqc_summary = join(OUT_DIR, 'multiqc', 'summary_files.txt'),
        multiqc_folder = join(OUT_DIR, 'multiqc')
    log:
        join(OUT_DIR, 'logs', 'MultiQC', 'multiQC.log')
    benchmark:
        join(OUT_DIR, 'benchmarks', 'MultiQC', 'multiQC.benchmark.tsv')
    message:
        """--- Running MultiQC for mapping"""
    conda:
        '../envs/vcf.yml'
    shell:
        """        
        ls -1 {params.fastqc_zip} >> {params.multiqc_summary}
        ls -1 {input.qc} >> {params.multiqc_summary}
        ls -1 {input.bamqc} | grep ":" | sed "s/://g" >> {params.multiqc_summary}
        multiqc -f -o {params.multiqc_folder} -d -dd 3 -l {params.multiqc_summary} > {log} 2>&1
        """