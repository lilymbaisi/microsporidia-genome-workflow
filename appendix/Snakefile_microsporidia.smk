configfile: "config.yaml"

rule all:
    input:
        expand("results/fastqc/{sample}_fastqc.html", sample=config["samples"]),
        expand("results/kraken2/{sample}_report.txt", sample=config["samples"]),
        expand("results/unicycler/{sample}/assembly.fasta", sample=config["samples"]),
        expand("results/busco/{sample}/short_summary.specific.microsporidia_odb10.txt", sample=config["samples"]),
        expand("results/augustus/{sample}.gff", sample=config["samples"])

rule fastqc:
    input:
        r1="data/{sample}_1.fq.gz",
        r2="data/{sample}_2.fq.gz"
    output:
        html1="results/fastqc/{sample}_1_fastqc.html",
        html2="results/fastqc/{sample}_2_fastqc.html"
    shell:
        "fastqc {input.r1} {input.r2} -o results/fastqc"

rule kraken2:
    input:
        r1="data/{sample}_1.fq.gz",
        r2="data/{sample}_2.fq.gz"
    output:
        report="results/kraken2/{sample}_report.txt",
        classified="results/kraken2/{sample}_classified#.fq",
        unclassified="results/kraken2/{sample}_unclassified#.fq"
    params:
        db=config["kraken2_db"]
    shell:
        "kraken2 --db {params.db} --paired --classified-out {output.classified} "
        "--unclassified-out {output.unclassified} --report {output.report} {input.r1} {input.r2}"

rule unicycler:
    input:
        r1="results/kraken2/{sample}_unclassified_1.fq",
        r2="results/kraken2/{sample}_unclassified_2.fq"
    output:
        assembly="results/unicycler/{sample}/assembly.fasta"
    shell:
        "unicycler -1 {input.r1} -2 {input.r2} -o results/unicycler/{wildcards.sample} --no_pilon --threads 32"

rule busco:
    input:
        fasta="results/unicycler/{sample}/assembly.fasta"
    output:
        summary="results/busco/{sample}/short_summary.specific.microsporidia_odb10.txt"
    params:
        lineage="microsporidia_odb10"
    shell:
        "busco -i {input.fasta} -o {wildcards.sample} -l {params.lineage} -m genome "
        "-f --out_path results/busco/{wildcards.sample} --download_path /mnt/lustre/bsp/DB/BUSCO/v5"

rule augustus:
    input:
        fasta="results/unicycler/{sample}/assembly.fasta"
    output:
        gff="results/augustus/{sample}.gff"
    params:
        species="microsporidia"
    shell:
        "augustus --species={params.species} {input.fasta} > {output.gff}"
