#!/usr/bin/env python3
"""
Yangtze River RNA Virus Pipeline - Snakemake Workflow
RNA virus bioinformatics pipeline for processing RNA NGS sequencing data
"""

# Configuration
configfile: "config/config_server/config_read.yaml"


# Helper to build optional CLI args for Ref_cluster_select.py
def build_ref_cluster_args(cfg):
    params = cfg.get("ref_cluster_select", {})
    args = []
    db_id_col = params.get("db_id_col")
    if db_id_col:
        args += ["--db-id-col", str(db_id_col)]
    sample_order = params.get("sample_order")
    if sample_order:
        args += ["--sample-order", str(sample_order)]
    sample_order_file = params.get("sample_order_file")
    if sample_order_file:
        args += ["--sample-order-file", str(sample_order_file)]
    if params.get("with_cluster_id"):
        args.append("--with-cluster-id")
    return " ".join(args)

REF_CLUSTER_ARGS = build_ref_cluster_args(config)


# Final output rule
rule all:
    input:
        "results/6_esviritu/databases/final_merged_database_only_assembly/virus_pathogen_database.all_metadata.tsv",
        "results/6_esviritu/databases/final_merged_database_only_assembly/virus_pathogen_database.fna",
        "results/6_esviritu/databases/final_merged_database_only_assembly/virus_pathogen_database.mmi",
        expand("results/6_esviritu/es/first_filter/{sample}/{sample}.detected_virus.assembly_summary.tsv", sample=config["samples"]),
        expand("results/6_esviritu/es/second_filter/{sample}/{sample}.detected_virus.assembly_summary.tsv", sample=config["samples"]),
        expand("results/6_esviritu/es/Merge/{sample}/{sample}.detected_virus.assembly_summary.merged.rpkmf.tsv", sample=config["samples"]),
        "results/6_esviritu/es/Merge/all_samples.detected_virus.assembly_summary.rpkmf.tsv",
        "results/6_esviritu/es/Correlation/all_samples.subspecies.min10.tsv",
        "results/6_esviritu/es/Correlation/all_samples.subspecies.spearman.all.tsv",
        "results/6_esviritu/es/Correlation/all_samples.subspecies.spearman.filtered.tsv",

# Rule 1: Quality control and adapter removal with fastp
# Note: fastp step already completed for current samples; outputs are in results/1_fastp/.
rule fastp_qc:
    input:
        r1 = "RNA_rawdata/{sample}/{sample}.R1.fq.gz",
        r2 = "RNA_rawdata/{sample}/{sample}.R2.fq.gz"
    output:
        r1_paired="results/1_fastp/{sample}/{sample}_1P.fq.gz",
        r2_paired="results/1_fastp/{sample}/{sample}_2P.fq.gz",
        r1_unpaired="results/1_fastp/{sample}/{sample}_U1.fq.gz",
        r2_unpaired="results/1_fastp/{sample}/{sample}_U2.fq.gz",
        html="results/1_fastp/{sample}/{sample}.fastp.html",
        json="results/1_fastp/{sample}/{sample}.fastp.json"     
    conda:
        "envs/fastp.yaml"
    resources:
        mem_mb_per_cpu = config["regular_memory"],  # MB
        runtime = config["fastp"]["runtime"],
        cpus_per_task = config["fastp"]["threads"],
        slurm_partition = config["regular_partition"],
        slurm_account = config["account"]
    log:
        out="log/1_fastp/{sample}.log",
        err="log/1_fastp/{sample}.err"
    shell:
        """
        mkdir -p results/1_fastp/{wildcards.sample}
        fastp --thread {threads} \
              --in1 {input.r1} --in2 {input.r2} \
              --out1 {output.r1_paired} --out2 {output.r2_paired} \
              --unpaired1 {output.r1_unpaired} --unpaired2 {output.r2_unpaired} \
              -h {output.html} -j {output.json} \
              --trim_poly_g --trim_poly_x \
              --qualified_quality_phred {config[fastp][quality_threshold]} \
              --length_required {config[fastp][length_required]} \
              --dont_overwrite > {log.out} 1> {log.err}  # Redirect stderr to log file
        """

# Rule 2: Remove rRNA sequences with ribodetector
rule ribodetector_rrna_removal:
    input:
        r1 = "results/1_fastp/{sample}/{sample}_1P.fq.gz",
        r2 = "results/1_fastp/{sample}/{sample}_2P.fq.gz"
    output:
        r1_nonrrna = "results/2_ribodetector/{sample}/{sample}_nonrrna.1.fq.gz",
        r2_nonrrna = "results/2_ribodetector/{sample}/{sample}_nonrrna.2.fq.gz"
    conda:
        "envs/ribodetector.yaml"
    resources:
        mem_mb_per_cpu = config["regular_memory"],  # MB
        runtime = config["ribodetector"]["runtime"],
        cpus_per_task = config["ribodetector"]["threads"],
        slurm_partition = config["regular_partition"],
        slurm_account = config["account"]
    log:
        out="log/2_ribodetector/{sample}.log",
        err="log/2_ribodetector/{sample}.err"
    shell:
        """
        mkdir -p results/2_ribodetector/{wildcards.sample}
        ribodetector_cpu -t {threads} \
                         -l {config[ribodetector][min_length]} \
                         -i {input.r1} {input.r2} \
                         -e rrna \
                         --chunk_size {config[ribodetector][chunk_size]} \
                         -o {output.r1_nonrrna} {output.r2_nonrrna} > {log.out} 1> {log.err}
        """
#Step 3: Assembly with SPAdes
rule spades_assembly:
    input:
        r1="results/2_ribodetector/{sample}/{sample}_nonrrna.1.fq.gz",
        r2="results/2_ribodetector/{sample}/{sample}_nonrrna.2.fq.gz"
    output:
        scaffolds="results/4_spades_result/{sample}/{sample}_no_correction/scaffolds.fasta"
    conda: "envs/spades.yaml"
    log:
        out="log/4_spades_assembly/{sample}_spades_assembly.log",
        err="log/4_spades_assembly/{sample}_spades_assembly.err"
    threads: config["spades"]["threads"]
    resources:
        slurm_partition=config["spades"]["spade_partition"],
        runtime=config["runtime"],
        mem_mb_per_cpu=int(config["spades"]["spade_memory"]),
        cpus_per_task=config["spades"]["threads"],
        slurm_account=config["slurm_account"]
    params:
        memory=lambda wildcards, resources: int(
            resources.mem_mb_per_cpu * resources.cpus_per_task / 1024
        ),
        k_list=config["spades"]["k_values"]
    shell:
        """
        mkdir -p results/4_spades_result/{wildcards.sample}
        mkdir -p log/4_spades_assembly
        spades.py --meta \
            -o results/4_spades_result/{wildcards.sample}/{wildcards.sample}_no_correction \
            -1 {input.r1} -2 {input.r2} \
            -t {threads} -m {params.memory} \
            -k {params.k_list} \
            --only-assembler \
            > {log.out} 2> {log.err}
        """

# Step 4: Rename and filter assemblies
rule rename_filter_assemblies:
    input:
        scaffolds="results/4_spades_result/{sample}/{sample}_no_correction/scaffolds.fasta"
    output:
        renamed=f"results/5_rename_assembly/rename_{config['seqkit']['min_length']}/{{sample}}_scaffolds_rename_{config['seqkit']['min_length']}.fasta"
    conda: "envs/seqkit.yaml"
    log:
        out="log/5_rename_filter/{sample}_rename_filter.log",
        err="log/5_rename_filter/{sample}_rename_filter.err"
    threads: config["seqkit"]["threads"]
    resources:
        slurm_partition=config["regular_partition"],
        runtime=config["runtime"],
        mem_mb_per_cpu=config["regular_memory"],
        cpus_per_task=config["seqkit"]["threads"],
        slurm_account=config["slurm_account"]
    params:
        min_length=config["seqkit"]["min_length"],
        nr_width=config["seqkit"]["nr_width"]
    shell:
        """
        mkdir -p results/5_rename_assembly/rename_{params.min_length}
        mkdir -p log/5_rename_filter


        # Filter sequences >= min_length and rename
        seqkit seq -m {params.min_length} {input.scaffolds} 2>> {log.err} | \
        seqkit replace -p .+ -r "{wildcards.sample}_{{nr}}" --nr-width {params.nr_width} \
        -o {output.renamed} \
        > {log.out} 2>> {log.err}
        """
# Step5.1 Merge with Esviritu database_just_merge
# Esviritu is a tool for identifying viral sequences from metagenomic data. however, After I try to mapping my filtered sequencing data to Esvirtu database, I found very few reads can be mapped to the database. So I decide to enlarge the database by merging established Esviritu database with my assembled contigs from metaspades. The minimum contig length is set to 200bp to ensure quality. which aligns with the Esviritu database construction criteria. for this step i just use seqkit to filter and cat all database together.

rule merge_esviritu_database_1:
    input:
        contigs = expand(
            f"results/5_rename_assembly/rename_{config['seqkit']['min_length']}/{{sample}}_scaffolds_rename_{config['seqkit']['min_length']}.fasta",
            sample=config["samples"],
        )
    output:
        merged_db = "results/6_esviritu/databases/database_merged_only_assembly/esviritu_merged_db_only_assembly.fasta"
    conda:
        "envs/seqkit.yaml"
    resources:
        mem_mb_per_cpu = config["regular_memory"],  # MB
        runtime = config["merge_esviritu_database_1"]["runtime"],  # minutes
        cpus_per_task = config["merge_esviritu_database_1"]["threads"],
        slurm_partition = config["regular_partition"],
        slurm_account = config["account"]
    params:
        min_length = config["merge_esviritu_database_1"]["min_length"]
    log:
        out="log/merge_esviritu_database/merge_esviritu_database_only_assembly_1.log",
        err="log/merge_esviritu_database/merge_esviritu_database_only_assembly_1.err"
    shell:
        """
        mkdir -p results/6_esviritu/databases/database_merged_only_assembly
        rm -f {output.merged_db}
        seqkit seq -m {params.min_length} {input.contigs} >> {output.merged_db}
        """


# # Merge with Esviritu database_step2/if have too much viral cluster
# # NOTE: Previously clustered with CheckV (BLAST + anicalc.py + aniclust.py).
# # Current strategy: mmseqs2 easy-cluster on nucleotide FASTA.
# # Step 15b: Cluster all samples using mmseqs2 method
# rule merge_esviritu_database_2:
#     input:
#         merged_db="results/6_esviritu/databases/database_merged/esviritu_merged_db.fasta"
#     output:
#         cluster_rep="results/6_esviritu/databases/cluster/all_samples_cluster_rep_seq.fasta",
#         cluster_all="results/6_esviritu/databases/cluster/all_samples_cluster_all_seqs.fasta",
#         cluster_tsv="results/6_esviritu/databases/cluster/all_samples_cluster.tsv"
#     conda: "envs/mmseqs2.yaml"
#     threads: config["mmseqs2"]["threads"]
#     log:
#         out="log/6_esviritu/mmseqs2_cluster.log",
#         err="log/6_esviritu/mmseqs2_cluster.err"
#     resources:
#         slurm_partition=config["regular_partition"],
#         runtime = config["mmseqs2"]["runtime"],
#         mem_mb_per_cpu = config["regular_memory"],
#         cpus_per_task = config["mmseqs2"]["threads"],
#         slurm_account = config["slurm_account"]
#     params:
#         min_seq_id=config["mmseqs2"]["min_seq_id"],
#         coverage=config["mmseqs2"]["coverage"],
#         tmp_dir="results/6_esviritu/databases/cluster/tmp",
#         prefix=config["mmseqs2"]["prefix"],
#         memory=config["regular_memory"]
#     shell:
#         """
#         mkdir -p results/6_esviritu/databases/cluster
#         mkdir -p log/6_esviritu
#         mkdir -p {params.tmp_dir}

#         # mmseqs2 easy-cluster on nucleotide FASTA
#         mmseqs easy-cluster {input.merged_db} \
#             {params.prefix} \
#             {params.tmp_dir} \
#             --min-seq-id {params.min_seq_id} \
#             -c {params.coverage} \
#             --threads {threads} > {log.out} 2> {log.err}
#         """

# Previous CheckV-based clustering strategy (kept for reference)
rule merge_esviritu_database_2_checkv:
    input:
        merged_db="results/6_esviritu/databases/database_merged_only_assembly/esviritu_merged_db_only_assembly.fasta"
    output:
        blast_db=directory("results/6_esviritu/databases/cluster_only_assembly/blast_db"),
        blast_results="results/6_esviritu/databases/cluster_only_assembly/all_samples_blast.tsv",
        ani_results="results/6_esviritu/databases/cluster_only_assembly/all_samples_ani.tsv",
        cluster_results="results/6_esviritu/databases/cluster_only_assembly/all_samples_cluster.tsv"
    conda: "envs/checkv.yaml"
    threads: config["checkv"]["threads"]
    log:
        out="log/6_esviritu/cluster_all.log",
        err="log/6_esviritu/cluster_all.err"
    resources:
        slurm_partition=config["checkv"]["slurm_partition"],
        runtime = config["checkv"]["runtime"],
        mem_mb_per_cpu = config["checkv"]["mem_mb_per_cpu"],
        cpus_per_task = config["checkv"]["threads"],
        slurm_account = config["slurm_account"]
    shell:
        """
        mkdir -p results/6_esviritu/databases/cluster_only_assembly
        mkdir -p log/6_esviritu

        # Create BLAST database
        makeblastdb -in {input.merged_db} \
            -dbtype nucl \
            -out results/6_esviritu/databases/cluster_only_assembly/all_samples_db \
            >> {log.out} 2>> {log.err}

        # Run BLAST all-vs-all
        blastn -query {input.merged_db} \
            -db results/6_esviritu/databases/cluster_only_assembly/all_samples_db \
            -outfmt '6 std qlen slen' \
            -max_target_seqs 1000 \
            -out {output.blast_results} \
            -num_threads {threads} \
            >> {log.out} 2>> {log.err}

        # Calculate ANI
        python {config[scripts][checkv_ani]} \
            -i {output.blast_results} \
            -o {output.ani_results} \
            >> {log.out} 2>> {log.err}

        # Cluster sequences
        python {config[scripts][checkv_clust]} \
            --fna {input.merged_db} \
            --ani {output.ani_results} \
            --out {output.cluster_results} \
            --min_ani {config[checkv][min_ani]} \
            --min_tcov {config[checkv][min_coverage]} \
            --min_qcov {config[checkv][min_qcov]} \
            >> {log.out} 2>> {log.err}

        # Create directory for blast database files
        mkdir -p {output.blast_db}
        mv results/6_esviritu/databases/cluster_only_assembly/all_samples_db.* {output.blast_db}/
        """


# Select representative sequences from each cluster
rule select_cluster_representatives:
    input:
        cluster="results/6_esviritu/databases/cluster_only_assembly/all_samples_cluster.tsv",
    output:
        reps="results/6_esviritu/databases/cluster_only_assembly/cluster_representatives.txt",
    conda:
        "envs/python.yaml"
    threads: 4
    resources:
        slurm_partition=config["regular_partition"],
        runtime=config["runtime"],
        mem_mb_per_cpu=config["regular_memory"],
        cpus_per_task=4,
        slurm_account=config["slurm_account"]
    log:
        out="log/6_esviritu/cluster_representatives.log",
        err="log/6_esviritu/cluster_representatives.err",
    shell:
        """
        mkdir -p log/6_esviritu
        python {config[scripts][ref_cluster_select_assembly_only]} \
            --cluster {input.cluster} \
            --output {output.reps} \
            {REF_CLUSTER_ARGS} \
            > {log.out} 2> {log.err}
        """

rule extract_final_virus_pathogen_database:
    input:
        reps="results/6_esviritu/databases/cluster_only_assembly/cluster_representatives.txt",
        fasta="results/6_esviritu/databases/database_merged_only_assembly/esviritu_merged_db_only_assembly.fasta",
    output:
        fasta="results/6_esviritu/databases/final_merged_database_only_assembly/virus_pathogen_database.fna",
    conda:
        "envs/seqkit.yaml"
    threads: config["seqkit"]["threads"]
    resources:
        slurm_partition=config["regular_partition"],
        runtime=config["runtime"],
        mem_mb_per_cpu=config["regular_memory"],
        cpus_per_task=config["seqkit"]["threads"],
        slurm_account=config["slurm_account"]
    log:
        out="log/6_esviritu/extract_final_db.log",
        err="log/6_esviritu/extract_final_db.err"
    shell:
        """
        mkdir -p results/6_esviritu/databases/final_merged_database_only_assembly
        mkdir -p log/6_esviritu
        seqkit grep -f {input.reps} {input.fasta} -o {output.fasta} \
            > {log.out} 2> {log.err}
        """

rule index_final_virus_pathogen_database:
    input:
        fasta="results/6_esviritu/databases/final_merged_database_only_assembly/virus_pathogen_database.fna",
    output:
        mmi="results/6_esviritu/databases/final_merged_database_only_assembly/virus_pathogen_database.mmi",
    conda:
        "envs/minimap2.yaml"
    threads: 4
    resources:
        slurm_partition=config["regular_partition"],
        runtime=config["runtime"],
        mem_mb_per_cpu=config["regular_memory"],
        cpus_per_task=4,
        slurm_account=config["slurm_account"]
    log:
        out="log/6_esviritu/minimap2_index.log",
        err="log/6_esviritu/minimap2_index.err"
    shell:
        """
        mkdir -p results/6_esviritu/databases/final_merged_database_only_assembly
        mkdir -p log/6_esviritu
        minimap2 -d {output.mmi} {input.fasta} \
            > {log.out} 2> {log.err}
        """

rule seqkit_length_from_merged_db:
    input:
        fasta="results/6_esviritu/databases/database_merged_only_assembly/esviritu_merged_db_only_assembly.fasta",
    output:
        length="results/6_esviritu/databases/final_merged_database_only_assembly/length.txt",
    conda:
        "envs/seqkit.yaml"
    threads: config["seqkit"]["threads"]
    resources:
        slurm_partition=config["regular_partition"],
        runtime=config["runtime"],
        mem_mb_per_cpu=config["regular_memory"],
        cpus_per_task=config["seqkit"]["threads"],
        slurm_account=config["slurm_account"]
    log:
        out="log/6_esviritu/seqkit_length.log",
        err="log/6_esviritu/seqkit_length.err"
    shell:
        """
        mkdir -p results/6_esviritu/databases/final_merged_database_only_assembly
        mkdir -p log/6_esviritu
        seqkit fx2tab -n -l {input.fasta} > {output.length} \
            2> {log.err}
        """

rule merge_esviritu_metadata:
    input:
        id_length="results/6_esviritu/databases/final_merged_database_only_assembly/length.txt",
        ids="results/6_esviritu/databases/cluster_only_assembly/cluster_representatives.txt",
        fasta="results/6_esviritu/databases/database_merged_only_assembly/esviritu_merged_db_only_assembly.fasta",
    output:
        merged="results/6_esviritu/databases/final_merged_database_only_assembly/virus_pathogen_database.all_metadata.tsv",
    conda:
        "envs/python.yaml"
    threads: 4
    resources:
        slurm_partition=config["regular_partition"],
        runtime=config["runtime"],
        mem_mb_per_cpu=config["regular_memory"],
        cpus_per_task=4,
        slurm_account=config["slurm_account"]
    log:
        out="log/6_esviritu/merge_metadata.log",
        err="log/6_esviritu/merge_metadata.err"
    shell:
        """
        mkdir -p results/6_esviritu/databases/final_merged_database_only_assembly
        mkdir -p log/6_esviritu
        python {config[scripts][merge_es_tsv_assembly_only]} \
            --id-length {input.id_length} \
            --ids {input.ids} \
            --fasta {input.fasta} \
            --output {output.merged} \
            > {log.out} 2> {log.err}
        """



rule esviritu_mapped_to_standard_database:
    input:
        r1="results/2_ribodetector/{sample}/{sample}_nonrrna.1.fq.gz",
        r2="results/2_ribodetector/{sample}/{sample}_nonrrna.2.fq.gz",
        db_dir=config["databases"]["esviritu_db_v3.2.4"],
    output:
        assembly_summary="results/6_esviritu/es/first_filter/{sample}/{sample}.detected_virus.assembly_summary.tsv",
        info="results/6_esviritu/es/first_filter/{sample}/{sample}.detected_virus.info.tsv",
        consensus="results/6_esviritu/es/first_filter/{sample}/{sample}_final_consensus.fasta",
        coverage="results/6_esviritu/es/first_filter/{sample}/{sample}.virus_coverage_windows.tsv",
        log_file="results/6_esviritu/es/first_filter/{sample}/{sample}_esviritu.log",
        params="results/6_esviritu/es/first_filter/{sample}/{sample}_esviritu.params.yaml",
        tagged_r1="results/6_esviritu/es/first_filter/{sample}/{sample}_nonrrna.1.tagged.fq.gz",
        tagged_r2="results/6_esviritu/es/first_filter/{sample}/{sample}_nonrrna.2.tagged.fq.gz",
        mapped_reads="results/6_esviritu/es/first_filter/{sample}/{sample}_temp/{sample}.reads.txt",
    conda:
        "envs/esviritu_map.yaml"
    threads: config["esviritu_map"]["threads"]
    resources:
        mem_mb_per_cpu=config["regular_memory"],
        runtime=config["esviritu_map"]["runtime"],
        cpus_per_task=config["esviritu_map"]["threads"],
        slurm_partition=config["regular_partition"],
        slurm_account=config["account"]
    params:
        temp_dir="results/6_esviritu/es/first_filter/{sample}/{sample}_temp",
        third_bam="results/6_esviritu/es/first_filter/{sample}/{sample}_temp/{sample}.third.filt.sorted.bam",
    log:
        out="log/6_esviritu/first_filter/{sample}.log",
        err="log/6_esviritu/first_filter/{sample}.err"
    shell:
        """
        mkdir -p results/6_esviritu/es/first_filter/{wildcards.sample}
        mkdir -p log/6_esviritu/first_filter
        zcat {input.r1} | gawk 'NR%4==1{{ $0=gensub(/^@([^ ]+)/, "@\\\\1_1", 1) }} {{ print }}' | gzip > {output.tagged_r1}
        zcat {input.r2} | gawk 'NR%4==1{{ $0=gensub(/^@([^ ]+)/, "@\\\\1_2", 1) }} {{ print }}' | gzip > {output.tagged_r2}
        python {config[scripts][modified_esviritu]} -r {output.tagged_r1} {output.tagged_r2} \
                 -s {wildcards.sample} \
                 -t {threads} \
                 -o results/6_esviritu/es/first_filter/{wildcards.sample} \
                 --db {input.db_dir} \
                 --keep True -q False \
                 > {log.out} 2>> {log.err}
        mkdir -p {params.temp_dir}
        samtools view {params.third_bam} | cut -f1 | sort -u > {output.mapped_reads}
        """

rule extract_unmapped_reads_esviritu:
    input:
        r1="results/6_esviritu/es/first_filter/{sample}/{sample}_nonrrna.1.tagged.fq.gz",
        r2="results/6_esviritu/es/first_filter/{sample}/{sample}_nonrrna.2.tagged.fq.gz",
        mapped_reads="results/6_esviritu/es/first_filter/{sample}/{sample}_temp/{sample}.reads.txt",
    output:
        unmapped_r1="results/6_esviritu/es/first_filter/{sample}/{sample}_nonrrna.1.unmapped.fq.gz",
        unmapped_r2="results/6_esviritu/es/first_filter/{sample}/{sample}_nonrrna.2.unmapped.fq.gz",
    conda:
        "envs/seqkit.yaml"
    threads: config["seqkit"]["threads"]
    resources:
        slurm_partition=config["regular_partition"],
        runtime=config["runtime"],
        mem_mb_per_cpu=config["regular_memory"],
        cpus_per_task=config["seqkit"]["threads"],
        slurm_account=config["slurm_account"]
    log:
        out="log/6_esviritu/unmapped_reads/{sample}.log",
        err="log/6_esviritu/unmapped_reads/{sample}.err"
    shell:
        """
        mkdir -p results/6_esviritu/es/first_filter/{wildcards.sample}
        mkdir -p log/6_esviritu/unmapped_reads
        tmp_dir=results/6_esviritu/es/first_filter/{wildcards.sample}/{wildcards.sample}_temp
        mkdir -p $tmp_dir
        seqkit seq -n {input.r1} > $tmp_dir/fastq_full_ids.txt
        seqkit seq -n {input.r2} >> $tmp_dir/fastq_full_ids.txt
        awk '
          NR==FNR {{
            core=$1
            full=$0
            map[core]=full
            next
          }}
          {{
            core=$1
            if (core in map) print map[core]
            else print core "\\tNOT_FOUND" > "/dev/stderr"
          }}
        ' $tmp_dir/fastq_full_ids.txt {input.mapped_reads} > $tmp_dir/ids.full.txt 2> {log.err}
        seqkit grep -v -n -f $tmp_dir/ids.full.txt {input.r1} -o {output.unmapped_r1} \
            >> {log.out} 2>> {log.err}
        seqkit grep -v -n -f $tmp_dir/ids.full.txt {input.r2} -o {output.unmapped_r2} \
            >> {log.out} 2>> {log.err}
        """

rule esviritu_mapped_to_assembly_database:
    input:
        r1="results/6_esviritu/es/first_filter/{sample}/{sample}_nonrrna.1.unmapped.fq.gz",
        r2="results/6_esviritu/es/first_filter/{sample}/{sample}_nonrrna.2.unmapped.fq.gz",
        db_dir="results/6_esviritu/databases/final_merged_database_only_assembly",
    output:
        assembly_summary="results/6_esviritu/es/second_filter/{sample}/{sample}.detected_virus.assembly_summary.tsv",
        info="results/6_esviritu/es/second_filter/{sample}/{sample}.detected_virus.info.tsv",
        consensus="results/6_esviritu/es/second_filter/{sample}/{sample}_final_consensus.fasta",
        coverage="results/6_esviritu/es/second_filter/{sample}/{sample}.virus_coverage_windows.tsv",
        log_file="results/6_esviritu/es/second_filter/{sample}/{sample}_esviritu.log",
        params="results/6_esviritu/es/second_filter/{sample}/{sample}_esviritu.params.yaml",
    conda:
        "envs/esviritu_map.yaml"
    threads: config["esviritu_map"]["threads"]
    resources:
        mem_mb_per_cpu=config["regular_memory"],
        runtime=config["esviritu_map"]["runtime"],
        cpus_per_task=config["esviritu_map"]["threads"],
        slurm_partition=config["regular_partition"],
        slurm_account=config["account"]
    params:
        extra=config["esviritu_map"].get("extra", "")
    log:
        out="log/6_esviritu/second_filter/{sample}.log",
        err="log/6_esviritu/second_filter/{sample}.err"
    shell:
        """
        mkdir -p results/6_esviritu/es/second_filter/{wildcards.sample}
        mkdir -p log/6_esviritu/second_filter
        python {config[scripts][modified_esviritu]} -r {input.r1} {input.r2} \
                 -s {wildcards.sample} \
                 -t {threads} \
                 -o results/6_esviritu/es/second_filter/{wildcards.sample} \
                 --db {input.db_dir} -q False \
                 {params.extra} \
                 > {log.out} 2>> {log.err}
        """

rule merge_esviritu_assembly_summary:
    input:
        first="results/6_esviritu/es/first_filter/{sample}/{sample}.detected_virus.assembly_summary.tsv",
        second="results/6_esviritu/es/second_filter/{sample}/{sample}.detected_virus.assembly_summary.tsv",
    output:
        first_copy="results/6_esviritu/es/Merge/{sample}/{sample}.detected_virus.assembly_summary.first.tsv",
        second_copy="results/6_esviritu/es/Merge/{sample}/{sample}.detected_virus.assembly_summary.second.tsv",
        merged="results/6_esviritu/es/Merge/{sample}/{sample}.detected_virus.assembly_summary.merged.tsv",
    resources:
        slurm_partition=config["regular_partition"],
        runtime=config["runtime"],
        mem_mb_per_cpu=config["regular_memory"],
        cpus_per_task=1,
        slurm_account=config["slurm_account"]
    log:
        out="log/6_esviritu/merge_info/{sample}.log",
        err="log/6_esviritu/merge_info/{sample}.err"
    shell:
        """
        mkdir -p results/6_esviritu/es/Merge/{wildcards.sample}
        mkdir -p log/6_esviritu/merge_info
        cp {input.first} {output.first_copy}
        cp {input.second} {output.second_copy}
        cat {input.first} {input.second} > {output.merged} 2> {log.err}
        """

rule rpkmf_esviritu_assembly_summary:
    input:
        merged="results/6_esviritu/es/Merge/{sample}/{sample}.detected_virus.assembly_summary.merged.tsv",
    output:
        rpkmf="results/6_esviritu/es/Merge/{sample}/{sample}.detected_virus.assembly_summary.merged.rpkmf.tsv",
    conda:
        "envs/python.yaml"
    resources:
        slurm_partition=config["regular_partition"],
        runtime=config["runtime"],
        mem_mb_per_cpu=config["regular_memory"],
        cpus_per_task=1,
        slurm_account=config["slurm_account"]
    log:
        out="log/6_esviritu/merge_info/{sample}.rpkmf.log",
        err="log/6_esviritu/merge_info/{sample}.rpkmf.err"
    shell:
        """
        python {config[scripts][esviritu_rpkmf]} --input {input.merged} --output {output.rpkmf} \
            > {log.out} 2> {log.err}
        """

rule merge_esviritu_rpkmf_all_samples:
    input:
        rpkmfs=expand(
            "results/6_esviritu/es/Merge/{sample}/{sample}.detected_virus.assembly_summary.merged.rpkmf.tsv",
            sample=config["samples"],
        ),
    output:
        merged="results/6_esviritu/es/Merge/all_samples.detected_virus.assembly_summary.rpkmf.tsv",
    conda:
        "envs/python.yaml"
    resources:
        slurm_partition=config["regular_partition"],
        runtime=config["runtime"],
        mem_mb_per_cpu=config["regular_memory"],
        cpus_per_task=1,
        slurm_account=config["slurm_account"]
    log:
        out="log/6_esviritu/merge_info/all_samples.rpkmf.log",
        err="log/6_esviritu/merge_info/all_samples.rpkmf.err"
    shell:
        """
        python {config[scripts][esviritu_merge_rpkmf]} --inputs {input.rpkmfs} --output {output.merged} \
            > {log.out} 2> {log.err}
        """

rule esviritu_subspecies_correlation:
    input:
        merged="results/6_esviritu/es/Merge/all_samples.detected_virus.assembly_summary.rpkmf.tsv"
    output:
        subspecies="results/6_esviritu/es/Correlation/all_samples.subspecies.min10.tsv",
        all_corr="results/6_esviritu/es/Correlation/all_samples.subspecies.spearman.all.tsv",
        filt_corr="results/6_esviritu/es/Correlation/all_samples.subspecies.spearman.filtered.tsv"
    threads:
        config.get("correlation", {}).get("threads", 8)
    conda:
        "envs/python.yaml"
    resources:
        slurm_partition=config["regular_partition"],
        runtime=config["runtime"],
        mem_mb_per_cpu=config["regular_memory"],
        cpus_per_task=1,
        slurm_account=config["slurm_account"]
    log:
        out="log/6_esviritu/correlation/subspecies_corr.log",
        err="log/6_esviritu/correlation/subspecies_corr.err"
    params:
        min_samples=config.get("correlation", {}).get("min_samples", 10),
        thresholds=config.get("correlation", {}).get("thresholds", "0.6,0.7,0.8,0.9"),
        p_threshold=config.get("correlation", {}).get("p_threshold", 0.05)
    shell:
        """
        mkdir -p results/6_esviritu/es/Correlation log/6_esviritu/correlation
        python scripts/esviritu_subspecies_correlation.py --input {input.merged} \
            --output-subspecies {output.subspecies} \
            --output-all {output.all_corr} \
            --output-filtered {output.filt_corr} \
            --min-samples {params.min_samples} \
            --thresholds {params.thresholds} \
            --p-threshold {params.p_threshold} \
            --threads {threads} \
            > {log.out} 2> {log.err}
        """
