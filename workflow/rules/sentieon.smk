__author__ = "Hanna Soderstrom"
__copyright__ = "Copyright 2023, Hanna Soderstrom"
__email__ = "hanna.soderstrom@gu.se"
__license__ = "GPL-3"


rule bwa_mem:
    input:
        reads=lambda wildcards: alignment_input(wildcards),
    output:
        bam=temp("sentieon/bwa_mem/{sample}_{flowcell}_{lane}_{barcode}_{type}.bam"),
    params:
        extra=config.get("sentieon", {}).get("extra", ""),
        reference=config.get("sentieon", {}).get("reference", ""),
        sentieon=config.get("sentieon", {}).get("sentieon", ""),
    log:
        "sentieon/bwa_mem/{sample}_{flowcell}_{lane}_{barcode}_{type}.output.log",
    benchmark:
        repeat(
            "sentieon/bwa_mem/{sample}_{flowcell}_{lane}_{barcode}_{type}.output.benchmark.tsv",
            config.get("sentieon", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("sentieon", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("sentieon", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("sentieon", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("sentieon", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("sentieon", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("sentieon", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("sentieon", {}).get("container", config["default_container"])
    conda:
        "../envs/sentieon.yaml"
    message:
        "{rule}: Align fastq files {input.reads} using Sentieon bwa mem against {params.reference}"
    shell:
        "{params.sentieon} bwa mem "
        "-M -R '@RG\\tID:{wildcards.sample}_{wildcards.type}\\tSM:{wildcards.sample}_{wildcards.type}\\tPL:ILLUMINA' "
        "-t {threads} {params.reference} {input.reads} "
        "| {params.sentieon} util sort -o {output.bam} "
        "-t {threads} "
        "--sam2bam "
        "-i - &> {log}"

rule dedup:
    input:
        lambda wildcards: [
            "sentieon/bwa_mem/{sample}_%s_%s_%s_{type}.bam" % (u.flowcell, u.lane, u.barcode)
            for u in get_units(units, wildcards)
        ],
    output:
        temp("sentieon/dedup/{sample}_{type}_DEDUP.bam"),
        "sentieon/dedup/{sample}_{type}_DEDUP.txt",
    params:
        extra=config.get("sentieon", {}).get("extra", ""),
        sentieon=config.get("sentieon", {}).get("sentieon", ""),
    log:
        "sentieon/dedup/{sample}_{type}.output.log",
    benchmark:
        repeat("sentieon/dedup/{sample}_{type}.output.benchmark.tsv", config.get("sentieon", {}).get("benchmark_repeats", 1))
    threads: config.get("sentieon", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("sentieon", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("sentieon", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("sentieon", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("sentieon", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("sentieon", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("sentieon", {}).get("container", config["default_container"])
    conda:
        "../envs/sentieon.yaml"
    message:
        "{rule}: Mark/remove duplicate reads in bam file {input} using Sentieon dedup algorithm"
    shell:
        "{params.sentieon} driver "
        "-t {threads} "
        "-i {input} "
        "--algo LocusCollector "
        "--fun score_info "
        "sentieon/dedup/{wildcards.sample}_{wildcards.type}_DEDUP_score.txt ; "
        "{params.sentieon} driver "
        "-t {threads} "
        "-i {input} "
        "--algo Dedup "
        "--rmdup "
        "--score_info sentieon/dedup/{wildcards.sample}_{wildcards.type}_DEDUP_score.txt "
        "--metrics sentieon/dedup/{wildcards.sample}_{wildcards.type}_DEDUP.txt "
        "sentieon/dedup/{wildcards.sample}_{wildcards.type}_DEDUP.bam &> {log}"

rule realigner:
    input:
        bam="sentieon/dedup/{sample}_{type}_DEDUP.bam",
    output:
        "sentieon/realign/{sample}_{type}_REALIGNED.bam",
    params:
        extra=config.get("sentieon", {}).get("extra", ""),
        reference=config.get("sentieon", {}).get("reference", ""),
        sentieon=config.get("sentieon", {}).get("sentieon", ""),
        mills=config.get("sentieon", {}).get("mills", ""),
    log:
        "sentieon/realign/{sample}_{type}.output.log",
    benchmark:
        repeat("sentieon/realign/{sample}_{type}.output.benchmark.tsv", config.get("sentieon", {}).get("benchmark_repeats", 1))
    threads: config.get("sentieon", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("sentieon", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("sentieon", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("sentieon", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("sentieon", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("sentieon", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("sentieon", {}).get("container", config["default_container"])
    conda:
        "../envs/sentieon.yaml"
    message:
        "{rule}: Indel realignment of bam file {input.bam} using Sentieon realigner"
    shell:
        "{params.sentieon} driver " 
        "-t {threads} " 
        "-r {params.reference} " 
        "-i {input.bam} " 
        "--algo Realigner " 
        "-k {params.mills} {output} &> {log}"

rule qualcal:
    input:
        bam="sentieon/realign/{sample}_{type}_REALIGNED.bam",
    output:
        "sentieon/qualcal/{sample}_{type}_RECAL_DATA.TABLE",
    params:
        extra=config.get("sentieon", {}).get("extra", ""),
        reference=config.get("sentieon", {}).get("reference", ""),
        sentieon=config.get("sentieon", {}).get("sentieon", ""),
        mills=config.get("sentieon", {}).get("mills", ""),
        dbsnp=config.get("sentieon", {}).get("dbsnp", ""),
    log:
        "sentieon/qualcal/{sample}_{type}.output.log",
    benchmark:
        repeat("sentieon/qualcal/{sample}_{type}.output.benchmark.tsv", config.get("sentieon", {}).get("benchmark_repeats", 1))
    threads: config.get("sentieon", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("sentieon", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("sentieon", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("sentieon", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("sentieon", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("sentieon", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("sentieon", {}).get("container", config["default_container"])
    conda:
        "../envs/sentieon.yaml"
    message:
        "{rule}: Calculate recalibration table of {input.bam} using Sentieon QualCal algorithm"
    shell:
        "{params.sentieon} driver "
        "-t {threads} "
        "-r {params.reference} "
        "-i {input.bam} "
        "--algo QualCal "
        "-k {params.mills} "
        "-k {params.dbsnp} {output} &> {log}"

rule dnascope:
    input:
        bam="sentieon/dedup/{sample}_{type}_DEDUP.bam",
    output:
        dnascope_vcf="sentieon/dnascope/{sample}_{type}_DNAscope.vcf",
        dnascope_idx="sentieon/dnascope/{sample}_{type}_DNAscope.vcf.idx",
    params:
        extra=config.get("sentieon", {}).get("extra", ""),
        reference=config.get("sentieon", {}).get("reference", ""),
        sentieon=config.get("sentieon", {}).get("sentieon", ""),
        callsettings=config.get("sentieon", {}).get("dnascope_settings", ""),
        model=config.get("sentieon", {}).get("dnascope_model", ""),
        dbsnp=config.get("sentieon", {}).get("dbsnp", ""),
    log:
        "sentieon/dnascope/{sample}_{type}.output.log",
    benchmark:
        repeat(
            "sentieon/dnascope/{sample}_{type}.output.benchmark.tsv",
            config.get("sentieon", {}).get("benchmark_repeats", 1)
        )
    threads: config.get("sentieon", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("sentieon", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("sentieon", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("sentieon", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("sentieon", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("sentieon", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("sentieon", {}).get("container", config["default_container"])
    conda:
        "../envs/sentieon.yaml"
    message:
        "{rule}: Call germline SNVs and structural variants in {input.bam} using Sentieon DNAScope"
    shell:
        "{params.sentieon} driver" 
        "-t {threads}" 
        "-r {params.reference} "
            "-i {input.bam}" 
            "--algo DNAscope" 
            "-d {params.dbsnp} "
            "--var_type snp,indel" 
            "--model {params.model} {params.callsettings} {output.dnascope_vcf}"

rule dnascope_modelfilter:
    input:
        vcf="sentieon/dnascope/{sample}_{type}_DNAscope.vcf",
        idx="sentieon/dnascope/{sample}_{type}_DNAscope.vcf.idx",
    output:
        vcf="sentieon/dnascope/{sample}_{type}_DNAscope_modelfiltered.vcf",
        idx="sentieon/dnascope/{sample}_{type}_DNAscope_modelfiltered.vcf.idx",
    params:
        extra=config.get("sentieon", {}).get("extra", ""),
        reference=config.get("sentieon", {}).get("reference", ""),
        sentieon=config.get("sentieon", {}).get("sentieon", ""),
        model=config.get("sentieon", {}).get("dnascope_model", ""),
    log:
        "sentieon/dnascope/{sample}_{type}_modelfiter.output.log",
    benchmark:
        repeat(
            "sentieon/dnascope/{sample}_{type}_modelfilter.output.benchmark.tsv",
            config.get("sentieon", {}).get("benchmark_repeats", 1)
        )
    threads: config.get("sentieon", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("sentieon", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("sentieon", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("sentieon", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("sentieon", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("sentieon", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("sentieon", {}).get("container", config["default_container"])
    conda:
        "../envs/sentieon.yaml"
    message:
        "{rule}: Modify the dnascope vcf {input.vcf} by adding the MLrejected filter to the variants using Sentieon DNAModelApply"
    shell:
        "{params.sentieon} driver -t {threads} "
        "-r {params.reference} "
        "--algo DNAModelApply "
        "--model {params.model} "
        "-v {input.vcf} {output.vcf} &> {log}"

rule tnscope:
    input:
        tumorbam="sentieon/realign/{sample}_T_REALIGNED.bam",
        normalbam="sentieon/realign/{sample}_N_REALIGNED.bam",
        tumortable="sentieon/qualcal/{sample}_T_RECAL_DATA.TABLE",
        normaltable="sentieon/qualcal/{sample}_N_RECAL_DATA.TABLE",
    output:
        tnscope="sentieon/tnscope/{sample}_TNscope_tn.vcf",
        tnscope_bam="sentieon/tnscope/{sample}_REALIGNED_realignedTNscope.bam",
    params:
        extra=config.get("sentieon", {}).get("extra", ""),
        reference=config.get("sentieon", {}).get("reference", ""),
        sentieon=config.get("sentieon", {}).get("sentieon", ""),
        callsettings=config.get("sentieon", {}).get("tnscope_settings", ""),
    log:
        "sentieon/tnscope/{sample}.output.log",
    benchmark:
        repeat("sentieon/tnscope/{sample}.output.benchmark.tsv", config.get("sentieon", {}).get("benchmark_repeats", 1))
    threads: config.get("tnscope", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("sentieon", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("sentieon", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("sentieon", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("tnscope", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("sentieon", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("sentieon", {}).get("container", config["default_container"])
    conda:
        "../envs/sentieon.yaml"
    message:
        "{rule}: Call SNVs and structural variants in {input.tumorbam} using matched normal {input.normalbam} using Sentieon TNScope"
    shell:
        "{params.sentieon} driver "
        "-t {threads} "
        "-r {params.reference} "
        "-i {input.tumorbam} "
        "-q {input.tumortable} "
        "-i {input.normalbam} "
        "-q {input.normaltable} "
        "--algo TNscope "
        "--tumor_sample {wildcards.sample}_T "
        "--normal_sample {wildcards.sample}_N "
        "--bam_output {output.tnscope_bam} "
        "{params.callsettings} {output.tnscope} &> {log}"

rule tnscope_modelfilter:
    input:
        tnscopevcf="sentieon/tnscope/{sample}_TNscope_tn.vcf",
        tnscopeidx="sentieon/tnscope/{sample}_TNscope_tn.vcf",
    output:
        vcf="sentieon/tnscope/{sample}_TNscope_tn_ML.vcf",
        idx="sentieon/tnscope/{sample}_TNscope_tn_ML.vcf.idx",
    params:
        extra=config.get("sentieon", {}).get("extra", ""),
        reference=config.get("sentieon", {}).get("reference", ""),
        sentieon=config.get("sentieon", {}).get("sentieon", ""),
        callsettings=config.get("sentieon", {}).get("tnscope_settings", ""),
        model=config.get("sentieon", {}).get("tnscope_model", ""),
    log:
        "sentieon/tnscope/{sample}.output.log",
    benchmark:
        repeat("sentieon/tnscope/{sample}.output.benchmark.tsv", config.get("sentieon", {}).get("benchmark_repeats", 1))
    threads: config.get("tnscope", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("sentieon", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("sentieon", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("sentieon", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("tnscope", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("sentieon", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("sentieon", {}).get("container", config["default_container"])
    conda:
        "../envs/sentieon.yaml"
    message:
        "{rule}: Apply machine learning model on the TNScope vcf {input.tnscopevcf} to help with variant filtration using Sentieon TNModelApply"
    shell:
        "{params.sentieon} driver "
        "-t {threads} "
        "-r {params.reference} "
        "--algo TNModelApply "
        "-m {params.model} "
        "-v {input.tnscopevcf} {output.vcf} &> {log}"
