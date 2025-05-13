fastqs = open('bams/StarAlign/samples-seq082024.txt').read().strip().split('\n')


rule all:
    input:
        expand("bams/StarAlign/CTRLs/{samples}/{samples}_Aligned.sortedByCoord.out.bam.bai", samples=fastqs)
        

rule STAR_align:
    input:
        fastq_R1="fastqs/RNAcap/{samples}_R1.fastq.gz",
        fastq_R2="fastqs/RNAcap/{samples}_R2.fastq.gz"
    output:
        bam="bams/StarAlign/CTRLs/{samples}/{samples}_Aligned.sortedByCoord.out.bam"
    log:
        "bams/StarAlign/tmp_output/{samples}.log"
    threads: 12
    resources:
        mem_gb=80
    params:
        genome_dir ="/mnt/gluster01/lgbm/Larissa.Souza/STAR/Genome-index-GRCh38/GenomeDir/",
        out_dir = "bams/StarAlign/CTRLs/{samples}/{samples}_",
        gtf_homo_sapiens = "/mnt/gluster01/lgbm/Larissa.Souza/anotacoes_homo_sapiens/Homo_sapiens.GRCh38.104.gtf",
        jobname = "Star_Align_{samples}"
    shell:
        'STAR --runThreadN 12 --genomeDir {params.genome_dir} --readFilesIn {input.fastq_R1} {input.fastq_R2} ' 
        '--readFilesCommand zcat --outFileNamePrefix {params.out_dir} --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif '
        '--outFilterMismatchNmax 3 --alignIntronMax 299999 --alignSJDBoverhangMin 6 --alignEndsType EndToEnd --chimSegmentMin 2 --sjdbGTFfile {params.gtf_homo_sapiens}'

rule samtools:
    input:
        bam="bams/StarAlign/CTRLs/{samples}/{samples}_Aligned.sortedByCoord.out.bam"
    output:
        bam_bai="bams/StarAlign/CTRLs/{samples}/{samples}_Aligned.sortedByCoord.out.bam.bai"
    log:
        "log/samtools/{samples}.log"
    threads: 8
    resources:
        mem_gb=10
    params:
        jobname = "samtools_{samples}"
    shell:
        '/conda/bin/samtools index -b {input} {output}'
        
