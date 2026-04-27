rule plink_ld:
    input:
        vg_vcf = join(OUT_DIR, 'calling', 'vg_annotated.vcf.gz'),
    output:
        plink_vcf = join(OUT_DIR, 'plink', 'vg_annotated_pruned.vcf')
    params:
        output_dir = directory(join(OUT_DIR, 'plink', 'ld')),
        output_prefix = join(OUT_DIR, 'plink', 'ld', 'vg_annotated_pruned'),
        ld_prune = join(OUT_DIR, 'plink', 'ld', 'vg_annotated_pruned.prune.in'),
        plink_vcf = join(OUT_DIR, 'plink', 'ld', 'vg_annotated_pruned.vcf')
    message:
        """--- Remove sites with linkage equilibrium."""
    conda:
        '../envs/plink.yml'
    shell:
        """
        mkdir -p {params.output_dir}
        plink --vcf {input.vg_vcf} --indep-pairwise 1000 1 0.2 --out {params.output_prefix} --double-id --allow-extra-chr
        plink --vcf {input.vg_vcf} --extract {params.ld_prune} --make-bed --out {params.output_prefix} --double-id --allow-extra-chr
        plink --bfile {params.output_prefix} --recode vcf --out {params.output_prefix} --double-id --allow-extra-chr
        mv {params.plink_vcf} {output.plink_vcf}
        rm -rf {params.output_dir}
        """

rule vcftools_stats:
    input:
        plink_vcf = join(OUT_DIR, 'plink', 'vg_annotated_pruned.vcf')
    output:
        join(OUT_DIR, 'plink_stats', 'vg_annotated_vcftools.log'),
    params:
        output_prefix = join(OUT_DIR, 'plink_stats', 'vg_annotated_vcftools'),
        window = 10000,
    message:
        """--- Calculate population statistics for pruned VCF."""
    conda:
        '../envs/plink.yml'
    shell:
        """
        vcftools --vcf {input.plink_vcf} --het --out {params.output_prefix}
        vcftools --vcf {input.plink_vcf} --hardy --out {params.output_prefix}
        vcftools --vcf {input.plink_vcf} --site-pi --out {params.output_prefix}
        vcftools --vcf {input.plink_vcf} --window-pi {params.window} --out {params.output_prefix}
        vcftools --vcf {input.plink_vcf} --fst-window-size {params.window} --out {params.output_prefix}
        """

rule vcftools_fst:
    input:
        plink_vcf = join(OUT_DIR, 'plink', 'vg_annotated_pruned.vcf'),
        fastq = expand(join(OUT_DIR, 'trimmed_reads', '{sample}.R1.fastq.gz'), sample = SAMPLES),
    output:
        join(OUT_DIR, 'plink_stats', 'vg_annotated_vcftools.weir.fst'),
    params:
        fastq_dir = join(OUT_DIR, 'trimmed_reads')
        fst_dir = temp(directory(join(OUT_DIR, 'plink_stats', 'fst'))),
        samples_list = temp(join(OUT_DIR, 'plink_stats', 'fst', 'samples.txt')),
        population_list = temp(join(OUT_DIR, 'plink_stats', 'fst', 'populations.txt')),
        output_prefix = join(OUT_DIR, 'plink_stats', 'vg_annotated_vcftools'),
    message:
        """--- Calculate Fst for pruned VCF."""
    conda:
        '../envs/plink.yml'
    shell:
        """
        mkdir {params.fst_dir}
        ls {params.fastq_dir} | sed '{{sub(".*/", "", $1)}} 1' | sed 's/\..*//g' | uniq > {params.samples_list}
        awk '{{prefix=$1; sub(/_.*$/, "", prefix); print $0, prefix}}' {params.samples_list} > {params.population_list}
        awk '{{print > $2 ".pop.txt"}}' {params.population_list}
        vcftools --vcf {input.plink_vcf} $(for f in {params.fst_dir}*.pop.txt; do printf -- "--weir-fst-pop %s " "$f"; done) --out {params.output_prefix}
        """