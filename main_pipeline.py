import os
import subprocess

def run_fastqc(reads_dir, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    subprocess.run(f"fastqc {reads_dir}/*.fq.gz -o {output_dir}", shell=True)
    subprocess.run(f"multiqc {output_dir}", shell=True)

def remove_host_reads(read1, read2, ref_genome, output_bam, threads=8):
    subprocess.run(f"bwa index {ref_genome}", shell=True)
    cmd = f"bwa mem -t {threads} {ref_genome} {read1} {read2} | " \
          f"samtools view -bS -f 4 - | samtools sort -o {output_bam}"
    subprocess.run(cmd, shell=True)

def run_kraken2(read1, read2, db, out_dir, sample):
    os.makedirs(out_dir, exist_ok=True)
    subprocess.run(f"kraken2 --db {db} --paired {read1} {read2} "
                   f"--classified-out {out_dir}/{sample}_classified#.fq "
                   f"--unclassified-out {out_dir}/{sample}_unclassified#.fq "
                   f"--report {out_dir}/{sample}_report.txt", shell=True)

def assemble_with_unicycler(read1, read2, out_dir):
    os.makedirs(out_dir, exist_ok=True)
    subprocess.run(f"unicycler -1 {read1} -2 {read2} -o {out_dir} --no_pilon --threads 32", shell=True)

def annotate_with_augustus(fasta, species, output_gff):
    subprocess.run(f"augustus --species={species} {fasta} > {output_gff}", shell=True)

def run_busco(fasta, lineage, output):
    subprocess.run(f"busco -i {fasta} -o {output} -l {lineage} -m genome "
                   f"--scaffold_composition -f --offline "
                   f"--download_path /mnt/lustre/bsp/DB/BUSCO/v5", shell=True)
