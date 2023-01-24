__author__ = "Hanna Soderstrom"
__copyright__ = "Copyright 2023, Hanna Soderstrom"
__email__ = "hanna.soderstrom@gu.se"
__license__ = "GPL-3"


rule bwa_mem:
    input:
        reads=lambda wildcards: alignment_input(wildcards),
    output:
        bam = "sentieon/bwa_mem/{sample}_{flowcell}_{lane}_{barcode}_{type}.bam", #FIXME Change to temp output
    params:
        extra=config.get("sentieon", {}).get("extra", ""),
        reference=config.get("sentieon", {}).get("reference", ""),
        sentieon=config.get("sentieon", {}).get("sentieon", ""),
    log:
        "sentieon/bwa_mem/{sample}_{flowcell}_{lane}_{barcode}_{type}.output.log",
    benchmark:
        repeat(
            "sentieon/bwa_mem/{sample}_{flowcell}_{lane}_{barcode}_{type}.output.benchmark.tsv",
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
        "{rule}: Align fastq files {input.reads} using Sentieon bwa mem against {params.reference}"
    shell:
        "{params.sentieon} bwa mem "
            "-M -R '@RG\\tID:{wildcards.sample}\\tSM:{wildcards.sample}\\tPL:ILLUMINA' "
            "-t {threads} {params.reference} {input.reads} "
        "| {params.sentieon} util sort -o {output.bam} -t {threads} --sam2bam -i -"

rule dedup:
    input:
        lambda wildcards: [
            "sentieon/bwa_mem/{sample}_%s_%s_%s_{type}.bam" % (u.flowcell, u.lane, u.barcode)
            for u in get_units(units, wildcards)
        ],
    output:
        "sentieon/dedup/{sample}_{type}_DEDUP.bam",
        "sentieon/dedup/{sample}_{type}_DEDUP.txt",
    params:
        extra=config.get("sentieon", {}).get("extra", ""),
        sentieon=config.get("sentieon", {}).get("sentieon", ""),
    log:
        "sentieon/dedup/{sample}_{type}.output.log",
    benchmark:
        repeat(
            "sentieon/dedup/{sample}_{type}.output.benchmark.tsv",
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
        "{rule}: Mark/remove duplicate reads in bam file {input} using Sentieon dedup algorithm"
    shell:
        "{params.sentieon} driver -t {threads} "
            "-i {input} "
            "--algo LocusCollector "
            "--fun score_info "
            "sentieon/dedup/{wildcards.sample}_{wildcards.type}_DEDUP_score.txt ;"
        "{params.sentieon} driver "
            "-t {threads} "
            "-i {input} "
            "--algo Dedup "
            "--rmdup "
            "--score_info sentieon/dedup/{wildcards.sample}_{wildcards.type}_DEDUP_score.txt "
            "--metrics sentieon/dedup/{wildcards.sample}_{wildcards.type}_DEDUP.txt "
            "sentieon/dedup/{wildcards.sample}_{wildcards.type}_DEDUP.bam"

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
        repeat(
            "sentieon/realign/{sample}_{type}.output.benchmark.tsv",
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
        "{rule}: Indel realignment of bam file {input.bam} using Sentieon realigner"
    shell:
        "{params.sentieon} driver -t {threads} -r {params.reference} -i {input.bam} --algo Realigner -k {params.mills} {output}"

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
        repeat(
            "sentieon/qualcal/{sample}_{type}.output.benchmark.tsv",
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
        "{rule}: Calculate recalibration table of {input.bam} using Sentieon QualCal algorithm"
    shell:
        "{params.sentieon} driver -t {threads} -r {params.reference} -i {input.bam} --algo QualCal -k {params.mills} -k {params.dbsnp} {output}"
