#Alinhamento para análise direcionada ao transcrito
bams = open('input-pipe.txt').read().strip().split('\n')
#bams='028_pool-lib-2053'

rule all:
    input:
        expand("Alinhamento-Star/{samples}/{samples}_Aligned.sortedByCoord.out.bam.bai", samples=bams)


rule bedtools:
    input:
        bam_ion="{samples}.bam"
    output:
        fastq="Fastqs/{samples}.fastq"
    log:
        "log/Fastqs/{samples}.log"
    threads: 12
    resources:
        mem_gb=80
    params:
        jobname = "Bedtools_{samples}"
    shell:
        '/conda/bin/bedtools bamtofastq -i {input} -fq {output}'

rule gzip:
    input:
        fastq="Fastqs/{samples}.fastq"
    output:
        fastq_gz="Fastqs/{samples}.fastq.gz"
    log:
        "log/Fastqs/{samples}_fq_gz.log"
    threads: 12
    resources:
        mem_gb=80
    params:
        jobname = "gzip_{samples}"
    shell:
        'gzip -c {input} > {output}'

rule STAR_align:
    input:
        fastq_gz="Fastqs/{samples}.fastq.gz"
    output:
        bam="Alinhamento-Star/{samples}/{samples}_Aligned.sortedByCoord.out.bam"
    log:
        "log/StarAlign/{samples}.log"
    threads: 12
    resources:
        mem_gb=80
    params:
        genome_dir ="/mnt/gluster01/lgbm/Larissa.Souza/STAR/Genome-index-GRCh38/GenomeDir/",
        out_dir = "Alinhamento-Star/{samples}/{samples}_",
        jobname = "Star_Align_{samples}",
        gtf_homo_sapiens = "/mnt/gluster01/lgbm/Larissa.Souza/anotacoes_homo_sapiens/Homo_sapiens.GRCh38.104.gtf"
    shell:
        'STAR --runThreadN 12 --genomeDir {params.genome_dir} --readFilesIn {input} --readFilesCommand zcat --outFileNamePrefix {params.out_dir} ' 
        '--outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --outFilterMismatchNmax 3 --alignIntronMax 299999 '
        '--alignSJDBoverhangMin 6 --alignEndsType EndToEnd --chimSegmentMin 2 --sjdbGTFfile {params.gtf_homo_sapiens}'

rule samtools:
    input:
        BAM_Star="Alinhamento-Star/{samples}/{samples}_Aligned.sortedByCoord.out.bam"
    output:
        BAM_BAI="Alinhamento-Star/{samples}/{samples}_Aligned.sortedByCoord.out.bam.bai"
    log:
        "log/samtools/{samples}.log"
    threads: 8
    resources:
        mem_gb=10
    params:
        jobname = "samtools_{samples}"
    shell:
        '/conda/bin/samtools index -b {input} {output}'