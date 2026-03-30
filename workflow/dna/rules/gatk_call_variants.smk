rule create_dict:
    input:
        dna = config['REFERENCE_FASTA']
    output:
        dict = os.path.splitext(config['REFERENCE_FASTA'])[0] + ".dict"
    log:
        join(OUT_DIR, 'logs', 'create_dict', 'dict.log')
    benchmark:
        join(OUT_DIR, 'benchmarks', 'create_dict', 'dict.benchmark.tsv')
    conda:
        '../envs/gatk.yml'
    shell:
        "gatk CreateSequenceDictionary -R {input.fasta} -O {output.dict} &> {log}"

rule GATK_allSites_calls:
    input:
        bam = join(OUT_DIR, 'alignment', '{sample}_alignment.bam'),
        dna = config['REFERENCE_FASTA'],
        dict = os.path.splitext(config['REFERENCE_FASTA'])[0] + ".dict"
    output:
        gatk_vcf = join(OUT_DIR, 'calling', 'gatk', '{sample}_gatk.vcf')
    log:
        join(OUT_DIR, 'logs', 'gatk', '{sample}.GATK_calls.log')
    benchmark:
        join(OUT_DIR, 'benchmarks', 'gatk', '{sample}.GATK_calls.benchmark.tsv')
    message:
        """--- Calling all sites with GATK for sample {wildcards.sample}."""
    threads:
        8
    resources:
        mem_mb=32000
    conda:
        '../envs/gatk.yml'
    shell:
        """
        gatk HaplotypeCaller -R {input.dna} -I {input.bam} -O {output.gatk_vcf} -ERC GVCF > {log} 2>&1
        """

rule GATK_gDBImport:
    input:
        gatk_vcf = expand(join(OUT_DIR, 'calling', 'gatk', '{sample}_gatk.vcf'), sample = SAMPLES)
    output:
        db_flag = join(OUT_DIR, 'calling', 'gatk', 'DBImport_complete')
    params:
        dna = config['REFERENCE_FASTA']
    log:
        join(OUT_DIR, 'logs', 'gatk', 'DBImport.log')
    benchmark:
        join(OUT_DIR, 'benchmarks', 'gatk', 'DBImport.benchmark.tsv')
    message:
        """--- Consolidating information from GVCF files across samples with GATK's GenomicsDBImport."""
    threads:
        8
    resources:
        mem_mb=32000
    conda:
        '../envs/gatk.yml'
    shell:
        """
        vcf_list=$(ls {input.gatk_vcf} | while read l; do echo " -V "$l; done)
        my_chr_list=$(grep ">" {params.dna} | sed "s/>//g" |  while read l; do echo " -L "$l; done)
        gatk GenomicsDBImport $vcf_list --genomicsdb-workspace-path allsamples_genomicsdb $my_chr_list > {log} 2>&1
        touch {output.db_flag}
        """

rule GATK_GenotypeGVCFs:
    input:
        db_flag = join(OUT_DIR, 'calling', 'gatk', 'DBImport_complete'),
        dna = config['REFERENCE_FASTA']
    output:
        gatk_vcf = join(OUT_DIR, 'calling', 'gatk.vcf.gz')
    log:
        join(OUT_DIR, 'logs', 'gatk', 'GenotypeGVCFs.log')
    benchmark:
        join(OUT_DIR, 'benchmarks', 'gatk', 'GenotypeGVCFs.benchmark.tsv')
    message:
        """--- Consolidating information from GVCF files across samples."""
    threads:
        24
    resources:
        mem_mb=256000
    conda:
        '../envs/gatk.yml'
    shell:
        """
        # Extract chromosome names from the reference genome
        chroms=$(grep ">" {input.dna} | sed 's/>//g')

        # Loop through each chromosome and run GenotypeGVCFs
        for chr in $chroms; do
            gatk GenotypeGVCFs \
                -R {input.dna} \
                -V gendb://allsamples_genomicsdb \
                -all-sites \
                -L $chr \
                -O {output.gatk_vcf}.$chr.vcf.gz \
                >> {log} 2>&1
        done

        # Concatenate all individual VCF files into one
        bcftools concat -o {output.gatk_vcf} -O z $(ls {output.gatk_vcf}.*.vcf.gz) >> {log} 2>&1
        """
