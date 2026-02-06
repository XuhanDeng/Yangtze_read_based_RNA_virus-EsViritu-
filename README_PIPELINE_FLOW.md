# Yangtze RNA Virus Pipeline (File Flow)

This is a simple, file‑centric overview of the current `Snakefile` workflow.

```
RAW READS
RNA_rawdata/{sample}/{sample}.R1/R2.fq.gz
   |
   v
fastp_qc
results/1_fastp/{sample}/{sample}_1P.fq.gz
results/1_fastp/{sample}/{sample}_2P.fq.gz
   |
   v
ribodetector_rrna_removal
results/2_ribodetector/{sample}/{sample}_nonrrna.1.fq.gz
results/2_ribodetector/{sample}/{sample}_nonrrna.2.fq.gz
   |
   +-------------------------> esviritu_map_final_database (reference DB)
   |                           DB dir: database/esviritu/v3.2.4
   |                           outputs per sample:
   |                           results/6_esviritu/es/first_filter/{sample}/
   |                             - {sample}.detected_virus.info.tsv
   |                             - {sample}.detected_virus.assembly_summary.tsv
   |                             - {sample}.virus_coverage_windows.tsv
   |                             - {sample}_final_consensus.fasta
   |                             - {sample}_esviritu.log + params.yaml
   |
   v
spades_assembly
results/4_spades_result/{sample}/{sample}_no_correction/scaffolds.fasta
   |
   v
rename_filter_assemblies
results/5_rename_assembly/rename_{minlen}/{sample}_scaffolds_rename_{minlen}.fasta
   |
   v
merge_esviritu_database_1 (all samples)
results/6_esviritu/databases/database_merged_only_assembly/
  esviritu_merged_db_only_assembly.fasta
   |
   v
merge_esviritu_database_2_checkv
results/6_esviritu/databases/cluster_only_assembly/
  - all_samples_blast.tsv
  - all_samples_ani.tsv
  - all_samples_cluster.tsv
  - blast_db/
   |
   v
select_cluster_representatives
results/6_esviritu/databases/cluster_only_assembly/
  cluster_representatives.txt
   |
   +--> extract_final_virus_pathogen_database
   |    results/6_esviritu/databases/final_merged_database_only_assembly/
   |      virus_pathogen_database.fna
   |
   +--> index_final_virus_pathogen_database
   |    results/6_esviritu/databases/final_merged_database_only_assembly/
   |      virus_pathogen_database.mmi
   |
   +--> seqkit_length_from_merged_db
   |    results/6_esviritu/databases/final_merged_database_only_assembly/
   |      length.txt
   |
   +--> merge_esviritu_metadata
        results/6_esviritu/databases/final_merged_database_only_assembly/
          virus_pathogen_database.all_metadata.tsv
```
