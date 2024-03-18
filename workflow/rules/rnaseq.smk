
sra_files, sra = glob_wildcards('rawdata/rnaseq/{sra_accession}/{sra_accession_}.sra')

KALLISTO_version = config["KALLISTO_version"] 
REF_BUILD = config["GENCODE_GENOME"]["REF_BUILD"]
REF_VERSION = config["GENCODE_GENOME"]["REF_VERSION"]




rule kallisto_quant:
    input:
        fastq1 = 'procdata/rnaseq/fastq/{sample}_1.fastq.gz',
        fastq2 = 'procdata/rnaseq/fastq/{sample}_2.fastq.gz',
        index = "references/{ref_build}_v{ref_version}/KALLISTO.v{KALLISTO_version}/kallisto-transcriptome.idx",
    output:
        files=multiext(
            "procdata/rnaseq/kallisto_v{KALLISTO_version}_{ref_build}.{ref_version}/{sample}/",
            "abundance.h5",
            "abundance.tsv",
            "run_info.json",
        ),
    log:
        'logs/rnaseq/kallisto_v{KALLISTO_version}_{ref_build}.{ref_version}/kallisto_quant/{sample}.log',
    envmodules: 
        f"kallisto/{KALLISTO_version}",
    threads: 8,
    shell:
        """
        kallisto quant \
            --threads {threads} \
            --index {input.index} \
            --output-dir $(dirname {output.files[0]}) \
            {input.fastq1} {input.fastq2} 2>&1 | tee {log}
        """

rule sra_to_fastq:
    input:
        sra = 'rawdata/rnaseq/{sra_accession}/{sra_accession}.sra',
    output:
        fastq1 = temp('procdata/rnaseq/fastq/{sra_accession}_1.fastq.gz'),
        fastq2 = temp('procdata/rnaseq/fastq/{sra_accession}_2.fastq.gz'),
    log:
        'logs/rnaseq/sra_to_fastq/{sra_accession}.log',
    envmodules:
        'sratoolkit/3.0.7',
        'pigz/2.6'
    threads: 8
    shell:
        """
        fasterq-dump \
            --split-3 \
            --progress \
            --details \
            --threads {threads} \
            --outdir procdata/rnaseq/fastq \
            --log-level info \
            {input.sra} 2>&1 | tee {log}

        echo "Compressing fastq files..." 2>&1 | tee -a {log}

        pigz --verbose \
            -p {threads} \
            -- procdata/rnaseq/fastq/{wildcards.sra_accession}_1.fastq \
            procdata/rnaseq/fastq/{wildcards.sra_accession}_2.fastq 2>&1 | tee -a {log}

        """