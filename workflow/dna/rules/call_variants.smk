rule prepare_pangenome:
    input:
        giraffe_gbz = join(GIRAFFE_FOLDER, '{prefix}.giraffe.gbz')
    output:
        full_xg = temp(join(OUT_DIR, 'calling', 'vg', '{prefix}.xg')),
        snarls = temp(join(OUT_DIR, 'calling', 'vg', '{prefix}.snarls'))
    threads:
        16
    message:
        """--- Annotate pangenome structure with VG. """
    conda:
        '../envs/vg.yml'
    shell:
        """
        vg convert -t {threads} -x {input.giraffe_gbz} > {output.full_xg}
        vg snarls -t {threads} {output.full_xg} > {output.snarls}
        """

rule vg_annotate:
    input: 
        vg_gam = join(OUT_DIR, 'trimmed_reads', '{sample}.gam'),
        full_xg = expand(join(OUT_DIR, 'calling', 'vg', '{prefix}.xg'), prefix = PREFIX),
        snarls = expand(join(OUT_DIR, 'calling', 'vg', '{prefix}.snarls'), prefix = PREFIX),
        gff = config['REFERENCE_GTF'],
    output:
        annotated_gam = temp(join(OUT_DIR, 'calling', 'vg', '{sample}.annotated.gam')),
    threads:
        16
    message:
        """--- Annotate '{wildcards.sample}' alignment with VG. """
    conda:
        '../envs/vg.yml'
    shell:
        """
        vg annotate -x {input.full_xg} -f {input.gff} -s {input.snarls} -a {input.vg_gam} -p -t {threads} > {output.annotated_gam}
        """

rule vg_variants:
    input:
        annotated_gam = join(OUT_DIR, 'calling', 'vg', '{sample}.annotated.gam'),
        full_xg = expand(join(OUT_DIR, 'calling', 'vg', '{prefix}.xg'), prefix = PREFIX),
        giraffe_gbz = expand(join(GIRAFFE_FOLDER, '{prefix}.giraffe.gbz'), prefix = PREFIX),
        snarls = expand(join(OUT_DIR, 'calling', 'vg', '{prefix}.snarls'), prefix = PREFIX)
    output:
        sample_vcf = temp(join(OUT_DIR, 'calling', 'vg', '{sample}_vg.vcf'))
    params:
        vg_pack = temp(join(OUT_DIR, 'calling', 'vg', '{sample}.pack')),
    threads:
        32
    message:
        """--- Call structural variants for sample {wildcards.sample} with VG. """
    conda:
        '../envs/vg.yml'
    shell:
        """
        vg pack -x {input.full_xg} -g {input.annotated_gam} -Q 5 -t {threads} -o {params.vg_pack}
        vg call -t {threads} -k {params.vg_pack} -r {input.snarls} -s {wildcards.sample} -a -O {input.giraffe_gbz} > {output.sample_vcf}
        """

rule vg_zip:
    input:
        sample_vcf = join(OUT_DIR, 'calling', 'vg', '{sample}_vg.vcf')
    output:
        vcf_zip = temp(join(OUT_DIR, 'calling', 'vg', '{sample}_vg.vcf.gz')),
        vcf_csi = temp(join(OUT_DIR, 'calling', 'vg', '{sample}_vg.vcf.gz.csi'))
    threads:
        32
    message:
        """--- Zip VG VCF for BCFTools. """
    conda:
        '../envs/bcftools.yml'
    shell:
        """
        bgzip {input.sample_vcf} -c > {output.vcf_zip}
        bcftools index {output.vcf_zip}
        """

rule vg_merge:
    input:
        vcf_zip = expand(join(OUT_DIR, 'calling', 'vg', '{sample}_vg.vcf.gz'), sample = SAMPLES),
        vcf_csi = expand(join(OUT_DIR, 'calling', 'vg', '{sample}_vg.vcf.gz.csi'), sample = SAMPLES),
    output:
        vg_vcf = join(OUT_DIR, 'calling', 'vg_annotated.vcf.gz')
    threads:
        16
    message:
        """--- Combining variant calls from VG. """
    conda:
        '../envs/bcftools.yml'
    shell:
        """
        bcftools merge -m all --force-samples -o {output.vg_vcf} --threads {threads} -O z {input.vcf_zip}
        bcftools +fill-tags {output.vg_vcf}  -- -t AF
        """

rule snpeffect:
    input:
        vg_vcf = join(OUT_DIR, 'calling', 'vg_annotated.vcf.gz'),
        snpeffect_script = config['SNPEFFECT_SCRIPT'],
        snpeffect_db = config['SNPEFFECT_DB']
    output:
        snpeffect_vcf = join(OUT_DIR, 'snpeffect', 'vg_annotated.snpEff.vcf.gz')
    params:
        snpeffect_vcf = join(OUT_DIR, 'snpeffect', 'vg_annotated.snpEff.vcf')
    threads:
        16
    message:
        """--- Predict variant effects from VG VCF with SnpEffect. """
    shell:
        """
        java -jar {input.snpeffect_script} -v {input.snpeffect_db} {input.vg_vcf} > {params.snpeffect_vcf}
        bgzip {params.snpeffect_vcf}
        """