Raw RNA-seq reads
├── RNA_rawdata/{sample}.R1.fq.gz
├── RNA_rawdata/{sample}.R2.fq.gz
│
├── fastp (QC + trimming)
│   ├── {sample}_1P.fq.gz
│   ├── {sample}_2P.fq.gz
│   ├── {sample}_U1.fq.gz
│   └── {sample}_U2.fq.gz
│
├── ribodetector (rRNA removal)
│   ├── {sample}_nonrrna.1.fq.gz
│   └── {sample}_nonrrna.2.fq.gz
│
├── ── Round 1: EsViritu (standard database) ──
│   │
│   ├── Tag read IDs
│   │   ├── append "_1" → {sample}_nonrrna.1.tagged.fq.gz
│   │   └── append "_2" → {sample}_nonrrna.2.tagged.fq.gz
│   │
│   ├── EsViritu mapping + filtering
│   │   ├── Reference DB:
│   │   │   └── esviritu_db_v3.2.4
│   │   │       ├── virus_pathogen_database.fna
│   │   │       └── virus_pathogen_database.mmi
│   │   │
│   │   ├── Outputs (Round 1 results)
│   │   │   ├── {sample}.detected_virus.info.tsv
│   │   │   ├── {sample}.detected_virus.assembly_summary.tsv
│   │   │   ├── {sample}_final_consensus.fasta
│   │   │   ├── {sample}.virus_coverage_windows.tsv
│   │   │   ├── {sample}_esviritu.log
│   │   │   └── {sample}_esviritu.params.yaml
│   │   │
│   │   └── Temporary BAM
│   │       └── {sample}.third.filt.sorted.bam
│   │           │
│   │           └── Extract mapped read names
│   │               └── samtools view | cut -f1 | sort -u
│   │                   │
│   │                   └── {sample}.reads.txt
│   │
│   ├── Reverse extraction (name-based)
│   │   ├── from tagged R1
│   │   │   └── seqkit grep -v
│   │   │       └── {sample}_nonrrna.1.unmapped.fq.gz
│   │   │
│   │   └── from tagged R2
│   │       └── seqkit grep -v
│   │           └── {sample}_nonrrna.2.unmapped.fq.gz
│
├── ── Round 2: EsViritu (assembly-derived database) ──
│   │
│   ├── Input reads
│   │   ├── {sample}_nonrrna.1.unmapped.fq.gz
│   │   └── {sample}_nonrrna.2.unmapped.fq.gz
│   │
│   ├── Reference DB
│   │   └── final_merged_database_only_assembly
│   │       ├── virus_pathogen_database.fna
│   │       ├── virus_pathogen_database.mmi
│   │       └── virus_pathogen_database.all_metadata.tsv
│   │
│   ├── EsViritu mapping + filtering
│   │
│   └── Outputs (Round 2 results)
│       ├── {sample}.detected_virus.info.tsv
│       ├── {sample}.detected_virus.assembly_summary.tsv
│       ├── {sample}_final_consensus.fasta
│       ├── {sample}.virus_coverage_windows.tsv
│       ├── {sample}_esviritu.log
│       └── {sample}_esviritu.params.yaml
│
└── Downstream analyses
    ├── Known virus profiles (Round 1)
    ├── Putative novel virus profiles (Round 2)
    └── Cross-sample correlation / ecology analysis