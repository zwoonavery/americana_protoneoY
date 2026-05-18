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

rule vcftools_heterozygosity:
    input:
        plink_vcf = join(OUT_DIR, 'plink', 'vg_annotated_pruned.vcf')
    output:
        join(OUT_DIR, 'plink_stats', 'vg_annotated_vcftools.het'),
    params:
        output_prefix = join(OUT_DIR, 'plink_stats', 'vg_annotated_vcftools'),
    message:
        """--- Calculate heterozygosity for pruned VCF."""
    conda:
        '../envs/plink.yml'
    shell:
        """
        vcftools --vcf {input.plink_vcf} --het --out {params.output_prefix}
        """

rule vcftools_hwe:
    input:
        plink_vcf = join(OUT_DIR, 'plink', 'vg_annotated_pruned.vcf')
    output:
        join(OUT_DIR, 'plink_stats', 'vg_annotated_vcftools.hwe'),
    params:
        output_prefix = join(OUT_DIR, 'plink_stats', 'vg_annotated_vcftools'),
    message:
        """--- Calculate Hardy-Weinberg Equilibrium for pruned VCF."""
    conda:
        '../envs/plink.yml'
    shell:
        """
        vcftools --vcf {input.plink_vcf} --hardy --out {params.output_prefix}
        """

rule vcftools_sites_pi:
    input:
        plink_vcf = join(OUT_DIR, 'plink', 'vg_annotated_pruned.vcf')
    output:
        join(OUT_DIR, 'plink_stats', 'vg_annotated_vcftools.sites.pi'),
    params:
        output_prefix = join(OUT_DIR, 'plink_stats', 'vg_annotated_vcftools'),
    message:
        """--- Calculate nucleotide diversity by site for pruned VCF."""
    conda:
        '../envs/plink.yml'
    shell:
        """
        vcftools --vcf {input.plink_vcf} --site-pi --out {params.output_prefix}
        """

rule vcftools_window_pi:
    input:
        plink_vcf = join(OUT_DIR, 'plink', 'vg_annotated_pruned.vcf')
    output:
        join(OUT_DIR, 'plink_stats', 'vg_annotated_vcftools.windowed.pi'),
    params:
        output_prefix = join(OUT_DIR, 'plink_stats', 'vg_annotated_vcftools'),
        pi_window = 10000,
    message:
        """--- Calculate nucleotide diversity by window for pruned VCF."""
    conda:
        '../envs/plink.yml'
    shell:
        """
        vcftools --vcf {input.plink_vcf} --window-pi {params.pi_window} --out {params.output_prefix}
        """

rule vcftools_pi_populations:
    input:
        fastq = expand(join(OUT_DIR, 'trimmed_reads', '{sample}.R1.fastq.gz'), sample = SAMPLES),
        plink_vcf = join(OUT_DIR, 'plink', 'vg_annotated_pruned.vcf')
    output:
        pi_population = join(OUT_DIR, 'plink_stats', 'vg_annotated_population_windowed.txt'),
    params:
        fastq_dir = join(OUT_DIR, 'trimmed_reads'),
        samples_list = temp(join(OUT_DIR, 'plink_stats', 'pi', 'samples.txt')),
        population_list = temp(join(OUT_DIR, 'plink_stats', 'pi', 'populations.txt')),
        pi_populations = temp(join(OUT_DIR, 'plink_stats', 'pi', '*.pop.txt')),
        output_prefix = join(OUT_DIR, 'plink_stats', 'vg_annotated_vcftools'),
        pi_window = 10000,
        pi_files = join(OUT_DIR, 'plink_stats', 'pi', '*.pop.txt.pi.txt')
        plink_combined = temp(join(OUT_DIR, 'plink_stats', 'pi', 'combined.txt')),
    message:
        """--- Calculate nucleotide diversity by window for pruned VCF."""
    conda:
        '../envs/plink.yml'
    shell:
        """
        ls {params.fastq_dir} | awk '{{sub(".*/", "", $1)}} 1' | sed 's/\..*//g' | uniq > {params.samples_list}
        awk '{{prefix=$1; sub(/_.*$/, "", prefix); print $0, prefix}}' {params.samples_list} > {params.population_list}
        awk '{{print > "{params.fst_dir}/" $2 ".pop.txt"}}' {params.population_list}
        for i in {params.pi_populations}; do
            vcftools --vcf {input.plink_vcf} --keep $i --window-pi {params.pi_window} --out $i
            sed '1d' $i.windowed.pi | awk 'BEGIN {{ FS=OFS="\t" }} {{print $0, "${{i}}"}}' > $i.pi.txt
        done
        cat {pi_files} >> {params.plink_combined}
        sed --expression '1i chrom\tbin_start\tbin_end\tn_variants\tn_monomorphic\tpi\tpop' {params.plink_combined} > {output.pi_population}
        """

rule vcftools_fst:
    input:
        plink_vcf = join(OUT_DIR, 'plink', 'vg_annotated_pruned.vcf'),
        fastq = expand(join(OUT_DIR, 'trimmed_reads', '{sample}.R1.fastq.gz'), sample = SAMPLES),
    output:
        plink_fst = join(OUT_DIR, 'plink_stats', 'vg_annotated_vcftools.weir.fst'),
    params:
        fastq_dir = join(OUT_DIR, 'trimmed_reads'),
        samples_list = temp(join(OUT_DIR, 'plink_stats', 'fst', 'samples.txt')),
        population_list = temp(join(OUT_DIR, 'plink_stats', 'fst', 'populations.txt')),
        fst_dir = temp(directory(join(OUT_DIR, 'plink_stats', 'fst'))),
        fst_populations = temp(join(OUT_DIR, 'plink_stats', 'fst', '*.pop.txt')),
        fst_window = 10000,
        fst_step = 1,
        plink_combined = temp(join(OUT_DIR, 'plink_stats', 'fst', 'combined.txt'))
    message:
        """--- Calculate Fst for pruned VCF."""
    conda:
        '../envs/plink.yml'
    shell:
        """
        ls {params.fastq_dir} | awk '{{sub(".*/", "", $1)}} 1' | sed 's/\..*//g' | uniq > {params.samples_list}
        awk '{{prefix=$1; sub(/_.*$/, "", prefix); print $0, prefix}}' {params.samples_list} > {params.population_list}
        awk '{{print > "{params.fst_dir}/" $2 ".pop.txt"}}' {params.population_list}
        populations=({params.fst_populations})
        for ((i = 0; i < ${{#populations[@]}}; i++)); do
            for ((j = i + 1; j < ${{#populations[@]}}; j++)); do
                vcftools --vcf {input.plink_vcf} --weir-fst-pop ${{populations[i]}} --weir-fst-pop ${{populations[j]}} --fst-window-size {params.fst_window} --fst-window-step {params.fst_step} --out ${{i}}_${{j}}
                sed '1d' ${{i}}_${{j}}.windowed.weir.fst | awk 'BEGIN {{ FS=OFS="\t" }} {{print $0, "${{populations[i]}}", ${{populations[j]}}}}' > ${{i}}_${{j}}.fst.txt
            done
        done
        cat *.fst.txt >> {params.plink_combined}
        sed --expression '1i chrom\tbin_start\tbin_end\tn_variants\tweighted_fst\tmean_fst\tpop1\tpop2' {params.plink_combined} > {output.plink_fst}
        """