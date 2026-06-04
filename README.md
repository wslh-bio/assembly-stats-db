# assembly-stats-db

This repository builds a per-species (per-taxid) statistics database of aggregated genome assembly statistics from the [NCBI RefSeq Genomes Assembly Summary File](https://ftp.ncbi.nlm.nih.gov/genomes/refseq/assembly_summary_refseq.txt)


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

This repository uses GitHub Actions to automatically generate a monthly update of the RefSeq assembly statistics database.

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

3. **Creates a Pull Request automatically**
    - A GitHub Actions workflow opens a PR containing the updated files

4. **Merge to main triggers release**
   - Once the PR is merged, a release workflow publishes a new versioned GitHub Release

---

## PR-Based Automation Model

This repository uses a pull request-based automation that ensures reproducibility and auditability of all database updates:

  - All automated updates are submitted as pull requests
  - The `main` branch is protected from direct pushes
  - Each update is reviewable and traceable
  - Releases are generated only after PR merge

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
│       ├── update-assembly-stats-db.yml
│       └── release-assembly-stats-db.yml
```

---

## Workflow Schedule

The database is updated automatically **once per month** using a scheduled GitHub Action:

```
schedule:
  - cron: "0 2 1 * *"
```  
This means the workflow will run at **02:00 UTC** on the **1st day of every month**

## Release behavior:

A GitHub Release is created automatically after each successful merge into `main`.
  - Trigger: push to `main`
  - Tag format: `vYYYYMMDD`
  - Assets: latest `NCBI_Assembly_Stats_YYYYMMDD.txt`
