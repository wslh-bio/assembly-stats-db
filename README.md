# assembly-stats-db

This repository is for a tool that builds a per-species (per-taxid) statistics database of aggregated genome assembly statistics from the [NCBI RefSeq Genomes Assembly Summary File](https://ftp.ncbi.nlm.nih.gov/genomes/refseq/assembly_summary_refseq.txt)


---

## Overview

The RefSeq assembly summary file contains per-assembly metadata such as genome size, GC content, CDS count, and taxonomy. This tool:

- Streams the RefSeq assembly summary file (local or URL)
- Aggregates values by **NCBI Taxonomy ID**
- Applies **IQR-based outlier filtering**
- Computes summary statistics per species
- Outputs a tab-delimited database

The resulting file follows the naming convention:

`NCBI_Assembly_Stats_YYYYMMDD.txt`

---

## Output Columns

| Column | Description |
|------|------------|
| Species | Scientific name |
| Min / Max / Median / Mean / StDev | Genome size statistics (Mb) |
| Assembly_count | Number of assemblies used |
| GC_Min / GC_Max / GC_Median / GC_Mean / GC_Stdev | GC content statistics (%) |
| GC_count | Number of assemblies with GC values |
| CDS_Min / CDS_Max / CDS_Median / CDS_Mean / CDS_Stdev | CDS count statistics |
| CDS_count | Number of assemblies with CDS values |
| Consensus_TAXID | NCBI taxonomy ID |

---

## Why Megabase Pairs (Mb)?

This tool reports genome sizes in **megabase pairs**, not base pairs.

Mb = base_pairs / 1,000,000

All genome size statistics are computed **after conversion**.

---

## IQR Filtering (Outlier Removal)

Assembly metadata can contain extreme values due to:

- Incomplete assemblies
- Contamination
- Mis-annotated records
- Mixed taxonomic assignments

To reduce distortion, this tool applies **Interquartile Range (IQR) filtering**:

- Values outside `Q1 − 1.5 × IQR` or `Q3 + 1.5 × IQR` are excluded
- Filtering is applied independently to:
  - Genome size
  - GC content
  - CDS count

This improves robustness of summary statistics while preserving biological variation.

Assemblies with fewer than four values are not filtered.

---

## Installation

Clone the repository:

```bash
git clone https://github.com/wslh-bio/assembly-stats-db.git
cd assembly-stats-db

```

Install dependencies:

`pip install numpy`

Python ≥ 3.8 is recommended


## Usage

## Output

The script generates the aggregated database named:

`NCBI_Assembly_Stats_YYYYMMDD.txt`