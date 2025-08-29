# Reference Data Files

This directory should contain the reference data files required for the TRACE pipeline.

## Required Files

### 1. hg38.fa.bed
RepeatMasker annotations for the human genome (hg38).

**Download and preparation:**
```bash
# Download from UCSC
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.out.gz

# Decompress
pigz -d hg38.fa.out.gz
# or
gunzip hg38.fa.out.gz

# Convert to BED format
python ../out2bed.py hg38.fa.out hg38.fa.bed
```

**File size:** ~306 MB (BED format)

### 2. Dfam-RepeatMasker.lib
Custom repeat library for RepeatMasker.

This file should be obtained from your specific analysis requirements or downloaded from:
- Dfam database: https://www.dfam.org/
- RepBase (requires subscription)
- Custom curated libraries

## Usage

When running TRACE.py, specify the paths to these files:

```bash
python ../TRACE.py input.vcf \
    --intersect data/hg38.fa.bed \
    --lib data/Dfam-RepeatMasker.lib \
    --threads 24
```

## Other Reference Genomes

For non-human genomes, obtain the corresponding RepeatMasker output files:

- **Mouse (mm10):** https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/
- **Zebrafish (danRer11):** https://hgdownload.soe.ucsc.edu/goldenPath/danRer11/bigZips/
- **Drosophila (dm6):** https://hgdownload.soe.ucsc.edu/goldenPath/dm6/bigZips/

Convert these to BED format using the same procedure as for hg38.

## Notes

- Large reference files are not included in the repository
- Users must download appropriate reference files for their organism
- Custom repeat libraries can be created using RepeatModeler or obtained from species-specific databases
