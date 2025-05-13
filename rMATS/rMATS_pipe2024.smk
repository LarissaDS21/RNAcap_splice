configfile: "rMATS_length.yaml"
BAMS = open('bams/StarAlign/Cases-seq082024.txt').read().strip().split('\n')
localrules: echo
Eventos=('A3SS','A5SS','MXE','RI','SE')

rule all:
    input:
        expand("results-2024/merged-rmats-tables/{samples}/{samples}_AS-filtrados.tsv", samples=BAMS)
        
rule echo:
    input:
        path_bam="bams/StarAlign/Casos/{samples}/{samples}_Aligned.sortedByCoord.out.bam"
    output: 
        b2_txt="arquivos-Casos-TXT/{samples}_b2.txt"
    shell:
        'echo {input} > {output}'

rule rMATS:
    input:
        b1_txt="arquivos-CTRLs-TXT/final_CTRLS2024.txt",
        b2_txt="arquivos-Casos-TXT/{samples}_b2.txt"
    output:
        A3SS_jc_txt="results-2024/{samples}/A3SS.MATS.JC.txt",
        A5SS_jc_txt="results-2024/{samples}/A5SS.MATS.JC.txt",
        MXE_jc_txt="results-2024/{samples}/MXE.MATS.JC.txt",
        RI_jc_txt="results-2024/{samples}/RI.MATS.JC.txt",
        SE_jc_txt="results-2024/{samples}/SE.MATS.JC.txt" 
    log:
        "results-2024/tmp_output/{samples}/{samples}_rmats.log"
    threads: 8
    resources:
        mem_gb=8
    params:
        out_dir = "results-2024/{samples}/",
        tmp_dir = "results-2024/tmp_output/{samples}/",
        length = lambda wc: config['length'][wc.samples],
        jobname = "rMATS_BAM_{samples}"
    shell:
        'rmats.py --b1 {input.b1_txt} --b2 {input.b2_txt} --gtf Homo_sapiens.GRCh38.104.gtf -t paired --variable-read-length --readLength {params.length} --novelSS '
        '--nthread {threads} --od {params.out_dir} --tmp {params.tmp_dir} > {log} 2>&1'

rule merged_tables:
    input:
        jc_txt=lambda wc:expand("results-2024/{samples}/{events}.MATS.JC.txt", samples=wc.samples, events=Eventos)
        # A3SS_jc_txt="results/{samples}/A3SS.MATS.JC.txt",
        # A5SS_jc_txt="results/{samples}/A5SS.MATS.JC.txt",
        # MXE_jc_txt="results/{samples}/MXE.MATS.JC.txt",
        # RI_jc_txt="results/{samples}/RI.MATS.JC.txt",
        # SE_jc_txt="results/{samples}/SE.MATS.JC.txt" 
    output:
        table_joined_nofilter="results-2024/merged-rmats-tables/{samples}/{samples}_noFilter.tsv",
        table_joined_filtered="results-2024/merged-rmats-tables/{samples}/{samples}_AS-filtrados.tsv"
    log:
        "results-2024/tmp_output/{samples}/{samples}_mergedtables.log"
    threads: 8
    resources:
        mem_gb=10
    params:
        jobname = "merged_tables_{samples}"
    conda:
        "r.yaml"
    notebook:
        "results-2024/template-merged-table.r.ipynb"