rule vg_vep:
    input:
        vcf_zip = join(OUT_DIR, 'calling', 'vg', '{sample}_vg.vcf.gz'),
        fasta_zip = config['FASTA_ZIP'],
        gtf_zip = config['GTF_ZIP'],
    output:
        vep_trim = join(OUT_DIR, 'vep_effect', '{sample}_vep_trim.txt')
    params:
        vep_out = temp(join(OUT_DIR, 'vep_effect', '{sample}_vep.txt')),
        vep_html = temp(join(OUT_DIR, 'vep_effect', '{sample}_vep.txt')),
        vep_summary = temp(join(OUT_DIR, 'vep_effect', '{sample}_vep.txt'))
    threads:
        16
    message:
        """--- Call variant effect for sample {wildcards.sample} with Ensembl's VEP. """
    conda:
        '../envs/vep.yml'
    shell:
        """
        vep -i {input.vcf_zip} --gff {input.gtf_zip} --fasta {input.fasta_zip} --fork {threads} --everything --force_overwrite -o {params.vep_out}
        awk '/>*/ {{ f = 1 }} f' {params.vep_out} | awk '{{print $0 "\t{wildcards.sample}"}}' > {output.vep_trim}
        """

rule vep_concatenate:
    input:
        vep_trim = expand(join(OUT_DIR, 'vep_effect', '{sample}_vep_trim.txt'), sample = SAMPLES)
    output:
        vep_cat = join(OUT_DIR, 'vep.tsv')
    message:
        """--- Concatenate Ensembl's VEP results. """
    shell:
        """
        cat {input.vep_trim} > {output.vep_cat}
        sed -i '1i #Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\tExisting_variation\tExtra\tSample' {output.vep_cat}
        """