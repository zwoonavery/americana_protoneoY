rule plink_ld:
    input:
        vg_vcf = join(OUT_DIR, 'calling', 'vg_annotated.vcf.gz'),
    output:
        plink_vcf = join(OUT_DIR, 'plink', 'vg_annotated_pruned.vcf')
    params:
        output_dir = directory(OUT_DIR, 'plink', 'ld'),
        output_prefix = join(OUT_DIR, 'plink', 'ld', 'vg_annotated_pruned'),
        ld_prune = join(OUT_DIR, 'plink', 'ld', 'vg_annotated_pruned.prune.in'),
        plink_bed = join(OUT_DIR, 'plink', 'ld', 'vg_annotated_pruned.bed'),
        plink_vcf = join(OUT_DIR, 'plink', 'ld', 'vg_annotated_pruned.vcf')
    message:
        """--- Remove sites with linkage equilibrium."""
    conda:
        '../envs/plink.yml'
    shell:
        """
        plink --vcf {input.vg_vcf} --indep 10000 1 0.2 --out {params.output_prefix} --double-id --allow-extra-chr
        plink --vcf (input.vg_vcf) --extract {params.ld_prune} --make-bed --out {params.output_prefix} --double-id --allow-extra-chr
        plink --bfile {params.plink_bed} --recode vcf --out {params.output_prefix} --double-id --allow-extra-chr
        mv {params.plink_vcf} {output_plink_vcf}
        rm -rf {params.output_dir}
        """

rule vcftools_stats:
    input:
        plink_vcf = join(OUT_DIR, 'plink', 'vg_annotated_pruned.vcf')
    output:
        join(OUT_DIR, 'plink', 'vg_annotated_vcftools.het'),
        join(OUT_DIR, 'plink', 'vg_annotated_vcftools.hwe'),
        join(OUT_DIR, 'plink', 'vg_annotated_vcftools.sites.pi'),
        join(OUT_DIR, 'plink', 'vg_annotated_vcftools.windowed.pi'),
        join(OUT_DIR, 'plink', 'vg_annotated_vcftools.weir.fst'),
    params:
        output_prefix = join(OUT_DIR, 'plink', 'vg_annotated_vcftools'),
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