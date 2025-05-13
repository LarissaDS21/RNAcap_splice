fastqs = open('/mnt/gluster01/lgbm/Larissa.Souza/PSI/DROP/RNAP_analyses/Bams/CasoseCTRLsinternos.txt').read().strip().split('\n')


rule all:
    input:
        expand("Casos/{samples}/{samples}_Aligned.sortedByCoord.out.bam.bai", samples=fastqs)
        

rule STAR_align:
    input:
        fastq_R1="/mnt/gluster01/lgbm/Larissa.Souza/PSI/rMATS/fastqs/RNAcap/{samples}_R1.fastq.gz",
        fastq_R2="/mnt/gluster01/lgbm/Larissa.Souza/PSI/rMATS/fastqs/RNAcap/{samples}_R2.fastq.gz"
    output:
        bam="Casos/{samples}/{samples}_Aligned.sortedByCoord.out.bam"
    log:
        "log/twopass_bams/{samples}.log"
    threads: 8
    resources:
        mem_gb=98
    params:
        genome_dir = "/mnt/gluster01/lgbm/Larissa.Souza/STAR/Genome-index-GRCh38-DROP/GenomeDir/",
        out_dir = "Casos/{samples}/{samples}_",
        jobname = "twopass_align_{samples}",
        gtf_homo_sapiens = "/mnt/gluster01/lgbm/Larissa.Souza/anotacoes_homo_sapiens/gencode.v46.primary_assembly.annotation.gtf"
    shell:
        'STAR --runThreadN 12 --genomeDir {params.genome_dir} --twopassMode Basic --readFilesIn {input.fastq_R1} {input.fastq_R2} ' 
        '--readFilesCommand zcat --outFileNamePrefix {params.out_dir} --outSAMtype BAM SortedByCoordinate '
        '--sjdbGTFfile {params.gtf_homo_sapiens}'

rule samtools_index:
    input:
        bam="Casos/{samples}/{samples}_Aligned.sortedByCoord.out.bam"
    output:
        bam_bai="Casos/{samples}/{samples}_Aligned.sortedByCoord.out.bam.bai"
    log:
        "log/samtools_index/{samples}.log"
    threads: 5
    resources:
        mem_gb=10
    params:
        jobname = "samtools_index_{samples}"
    shell:
        '/conda/bin/samtools index {input} {output}'