import os

juncfiles = open('wildcards_all-pipe.txt').read().strip().split('\n')
#conda_env:/home/25313@corp.accamargo.org.br/leafcutter
#export PATH="/home/25313@corp.accamargo.org.br/leafcutter-test/leafcutter/scripts:$PATH"
#export PATH="/conda/bin:$PATH"

rule all:
    input:
        expand("outlier_splice/out_24/{samples}/ctrlsVScase_clusterPvals.txt", samples=juncfiles)

rule junc_txt:
    input:
        patient_junc="Junctions/all-final/{samples}_Aligned.sortedByCoord.out.bam.junc"
    output: 
        patient_junc_txt="{samples}_juncfile.txt"
    log:
        "log/junc_txt/{samples}_juncfile.log"
    threads: 8
    resources:
        mem_gb=10
    params:
        jobname = "junc_txt_{samples}"
    shell:
        "echo {input.patient_junc} > {output}"
                
rule cat_jun_txt:
    input:
        patient_junc_txt="{samples}_juncfile.txt",
        ctrls_junc_txt="ctrls_srr-final.txt"
    output: 
        final_junc_txt="{samples}_final_juncfile.txt"
    log:
        "log/cat_jun_txt/{samples}_cat_jun_txt.log"
    threads: 8
    resources:
        mem_gb=10
    params:
        jobname = "cat_jun_txt_{samples}"
    shell:
        "cat {input.patient_junc_txt} {input.ctrls_junc_txt} > {output}"

rule intron_clustering:
    input:
        final_junc_txt="{samples}_final_juncfile.txt"
    output: 
        perind_numers_counts="intron_clustering/out_24/{samples}/ctrlsVScase_perind_numers.counts.gz"
    log:
        "log/intron_clustering/{samples}.log"
    threads: 8
    resources:
        mem_gb=10
    params:
        jobname = "intron_clustering_{samples}",
        outdir = lambda wc, output: os.path.dirname(output.perind_numers_counts)
    shell:
        "python leafcutter/clustering/leafcutter_cluster.py -j {input} -m 50 -o ctrlsVScase -l 500000 -r {params.outdir}"

rule Outlier_intron:
    input:
        perind_numers_counts="intron_clustering/out_24/{samples}/ctrlsVScase_perind_numers.counts.gz"
    output: 
        outlier_pVals="outlier_splice/out_24/{samples}/ctrlsVScase_clusterPvals.txt"
    log:
        "log/Outlier_intron/{samples}.log"
    threads: 8
    resources:
        mem_gb=10
    params:
        jobname = "Outlier_intron_{samples}",
        prefix = "outlier_splice/out_24/{samples}/ctrlsVScase"
    shell:
        "leafcutter/scripts/leafcutterMD.R --num_threads {threads} {input} -o {params.prefix}"