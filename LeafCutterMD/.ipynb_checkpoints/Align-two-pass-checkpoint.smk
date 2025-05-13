fastqs = open('twopass_bams/samples-seq082024.txt').read().strip().split('\n')


rule all:
    input:
        expand("twopass_bams/Cases/{samples}/{samples}_Aligned.sortedByCoord.out.bam.bai", samples=fastqs)
        

rule STAR_align:
    input:
        fastq_R1="/mnt/gluster01/lgbm/Larissa.Souza/PSI/rMATS/fastqs/RNAcap/{samples}_R1.fastq.gz",
        fastq_R2="/mnt/gluster01/lgbm/Larissa.Souza/PSI/rMATS/fastqs/RNAcap/{samples}_R2.fastq.gz"
    output:
        bam="twopass_bams/Cases/{samples}/{samples}_Aligned.sortedByCoord.out.bam"
    log:
        "log/twopass_bams/{samples}.log"
    threads: 8
    resources:
        mem_gb=98
    params:
        genome_dir = "/mnt/gluster01/lgbm/Larissa.Souza/STAR/Genome-index-GRCh38/GenomeDir/",
        out_dir = "twopass_bams/Cases/{samples}/{samples}_",
        jobname = "twopass_align_{samples}",
        gtf_homo_sapiens = "/mnt/gluster01/lgbm/Larissa.Souza/anotacoes_homo_sapiens/Homo_sapiens.GRCh38.104.gtf"
    shell:
        'STAR --runThreadN 12 --genomeDir {params.genome_dir} --twopassMode Basic --readFilesIn {input.fastq_R1} {input.fastq_R2} ' 
        '--readFilesCommand zcat --outFileNamePrefix {params.out_dir} --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif '
        '--sjdbGTFfile {params.gtf_homo_sapiens}'

rule samtools:
    input:
        bam="twopass_bams/Cases/{samples}/{samples}_Aligned.sortedByCoord.out.bam"
    output:
        bam_bai="twopass_bams/Cases/{samples}/{samples}_Aligned.sortedByCoord.out.bam.bai"
    log:
        "log/samtools/{samples}.log"
    threads: 8
    resources:
        mem_gb=10
    params:
        jobname = "samtools_{samples}"
    shell:
        '/conda/bin/samtools index -b {input} {output}'