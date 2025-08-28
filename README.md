# TRACE / tepipe

### A one-click user-friendly pipeline to annotate structural variants with transposable element information
#### Hackathon team: Michal Izydorczyk, Thomas Garcia, John Adedeji, Wai Yee (Nicola) Wong, Julian Chiu,  Quang Tran

TRACE is a pipeline designed to annotate structural variants (SV) with transposable element (TE) information and generate both an annotated VCF and a human-readable summary report.

## TRACE Workflow
<img width="6480" height="2016" alt="image" src="https://github.com/user-attachments/assets/dd4f9e96-455c-4add-a9e9-91d00950f663" />

## Workflow method

A reference sample BAM file was obtained from NCBI at [https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/Ultralong_OxfordNanopore/guppy-V3.2.4_2020-01-22/. ](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385[…]2020-01-22/HG002_GRCh38_ONT-UL_GIAB_20200122.phased.bam).

Reference FASTA: GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta

## Installation

## Inputs
- BAM/CRAM files
- Reference genome (FASTA)
- Optional: pre-existing SV VCF

## Outputs
- Annotated VCF with TE fields
- Summary report (HTML/Markdown/PDF)

## Pipeline summary
### Read alignment

This step is optional but can be enabled by using `--align` (?) and providing a FASTQ file

### SV calling on ONT data

Default structural variant caller is Sniffle2(?), however different callers can also be used by enabling  `--snv-caller`:
  - `sniffle`
  - `cutesv`
  - `svim` ?

### Repetitive element annotation

### RM parsing

### SV-Repeat Merging + Reporting

## Usage

DRAFT

`tepipe run --bam sample.bam --ref hg38.fa --input-vcf-type vcf-type --sv-caller [sniffle/cutesv]`

or

```
trace run . \
  --bam sample.bam \
  --ref hg38.fa \
  --input-vcf-type vcf-type \
  --sv-caller <sniffle/cutesv>
```

## Tools in Dockerfile

[Minimap2 (version used?)](https://minimap2.com/) is a versatile sequence alignment program that aligns DNA or mRNA sequences against a large reference database.

[Sniffles2 (Version used?)](https://github.com/fritzsedlazeck/Sniffles) is a fast structural variant caller for long-read sequencing.

[cuteSV (Version used?)](https://github.com/tjiangHIT/cuteSV) is a long-read based human genomic structural variation detection tool.

[svim (Version used)](https://github.com/eldariont/svim) is a structural variant caller for third-generation sequencing reads.

[other SV caller (Version used?)](url) temp.

~[bcftools](https://samtools.github.io/bcftools/) is a set of utilities that manipulate variant calls in the Variant Call Format (VCF) and its binary counterpart BCF.~

[Dfam RepeatMasker v4.1.8](https://www.dfam.org/releases/current/families/Dfam-RepeatMasker.lib.gz) is a program that screens DNA sequences for interspersed repeats and low complexity DNA sequences.

[One code to find them all v1.0](https://doua.prabi.fr/software/one-code-to-find-them-all) is a set of perl scripts to extract useful information from RepeatMasker about transposable elements, retrieve their sequences and get some quantitative information.

[bedtools (Version used?)](https://bedtools.readthedocs.io/en/latest/) is a set of utilities for a wide-range of genomics analysis tasks.

[snakemake/nextflow]  is a workflow system tool for creating reproducible, portable and scalable data analyses.

[Docker](https://www.docker.com/) is a container service provider that enhances reproducibility. (Or Conda)
