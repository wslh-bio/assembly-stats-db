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

## Megabase Pairs (Mb)

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

`python3 bin/calculate_assembly_stats.py -d assets/summary/assembly_summary_refseq_YYYYMMDD.txt.gz`

## Output

The script generates the aggregated database named:

`NCBI_Assembly_Stats_YYYYMMDD.txt`


---

## Automated Monthly Database Updates (GitHub Actions)

This repository includes a GitHub Actions workflow that automatically updates the RefSeq assembly statistics database on a monthly basis.

The workflow performs the following steps:

1. **Downloads the latest NCBI RefSeq assembly summary file**
   - Source:
     ```
     https://ftp.ncbi.nlm.nih.gov/genomes/refseq/assembly_summary_refseq.txt
     ```
   - The file is compressed and stored locally as:

     ```
     assets/summary/assembly_summary_refseq_YYMMDD.txt.gz
     ```

2. **Generates an aggregated assembly statistics database**
   - Runs:

     ```bash
     python3 bin/calculate_assembly_stats.py -d assets/summary/assembly_summary_refseq_YYMMDD.txt.gz
     ```

   - Produces:

     ```
     assets/database/NCBI_Assembly_Stats_YYMMDD.txt
     ```

3. **Overwrites previous databases**
   - Only one summary file and one generated database are retained at any time.
   - Each monthly update replaces the previous version.

4. **Commits the updated files back to the repository**
   - Commits are created automatically by the GitHub Actions bot.

---

### Repository Layout After Automation

```
assembly-stats-db
├── assets
│   ├── summary
│   │   └── assembly_summary_refseq_YYMMDD.txt.gz
│   └── database
│       └── NCBI_Assembly_Stats_YYMMDD.txt
├── bin
│   └── calculate_assembly_stats.py
├── .github
│   └── workflows
│       └── update-assembly-stats.yml
```

---

## Workflow Schedule

The database is updated automatically **once per month** using a scheduled GitHub Action:

```
schedule:
  - cron: "0 2 1 * *"
  
This means the workflow will run at 02:00 UTC on the 1st day of every month

When the workflow runs, a new release tagged with the date in YYMMDD format will be generated.