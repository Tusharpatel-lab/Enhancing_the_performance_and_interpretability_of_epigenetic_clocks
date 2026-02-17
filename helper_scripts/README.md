Description: This repository contains bash scripts for downloading and processing ATAC-seq data, performing footprinting analysis via TOBIAS, and annotating transcription factor binding sites (TFBS). The overall workflow is designed utilizing the latest bashbone release (More: https://github.com/Hoffmann-Lab/bashbone) on HPC Linux systems.

Directory Contents:

    > TFBS_annotation.sh: Annotates whole genome transcription factor binding sites by intersecting JASPAR motifs with Cistrome ChIP-seq peaks. Input is combined ChIP peak file, and other standard genomic files. 

    > getting_BAMS.sh: Downloads and prepares BAM files required for downstream footprinting analysis. Input is the list of SRR IDs (as srr.list), available on SRA FTP server.

    > Footprinting.sh: Performs ATAC-seq footprinting analysis to identify putative transcription factor binding events. Input is BAM files generated from the last step, and other standard genomic files. 


Tools required (versions as implied by bashbone, unless specified): 

| Script               | Tool                           | Purpose                                                                       |
| -------------------- | ------------------------------ | ---------------------------------------------------------                     |
| `TFBS_annotation.sh` | **bigBedToBed** (09-05-2019)   | Convert from bigBed to ascii bed format (UCSC)                                |
|                      | **bedtools**                   | Overlap footprints with TF binding sites                                      |
|                      | **awk / sed / coreutils**      | Formatting and filtering annotation outputs                                   |
|                      | **bashbone**                   | Parallel execution and job orchestration                                      |
| `getting_BAMS.sh`    | **wget**                       | Download FASTQ / sequencing data                                              |
|                      | **pigz**                       | Parallel decompression of FASTQ files                                         |
|                      | **cutadapt**                   | Adapter trimming and read preprocessing                                       |
|                      | **bwa-mem2**                   | Align ATAC-seq reads to reference genome                                      |
|                      | **samtools**                   | BAM conversion, sorting, indexing, duplicate removal                          |
|                      | **awk / sed / coreutils**      | Formatting and filtering annotation outputs                                   |
|                      | **bashbone**                   | Command construction, parallelization, HPC job submission                     |
| `Footprinting.sh`    | **samtools**                   | BAM filtering and indexing                                                    |
|                      | **deepTools** (v3.5.5)         | Genomic interval operations (Bigwig operations, and computing signals)        |
|                      | **TOBIAS** (v0.17.3)           | Transcription factor Occupancy prediction By Investigation of ATAC-seq Signal |
|                      | **awk / sed / coreutils**      | Text and file processing                                                      |
|                      | **bashbone**                   | Workflow control and parallel execution                                       |

NOTE:Paths to reference genomes, indices, and software installations must be adjusted to your local environment.