#!/usr/bin/env python3
import argparse
import logging
import urllib.request
import numpy as np
import statistics
from datetime import datetime
import gzip
from collections import defaultdict


logging.basicConfig(level = logging.INFO, format = '%(levelname)s : %(message)s')

timestamp = datetime.now().strftime("%Y%m%d")

def parse_args(args=None):
    Description='Generate a database with assembly statistics from RefSeq assembly summary file. '

    parser = argparse.ArgumentParser(description=Description)
    parser.add_argument('-d', '--path_database',
        metavar='path_to_database_file', 
        type=str,
        help='Path or URL of the Refseq assembly summary file', 
        required=True
        )
    parser.add_argument('-V', '--version',
        action='store_true', 
        help='Print version and exit'
        )
    return parser.parse_args(args)


def open_gzip(path):
    """
    Open a gzipped file from a local path or URL.
    """
    if path.startswith(("http://", "https://")):
        return gzip.open(urllib.request.urlopen(path), "rt")
    return gzip.open(path, "rt")


def is_float(value):
    """
    Helper function to handle strings in columns of refseq assembly file.
    Returns True if the value can be converted to a float.
    """
    try:
        float(value)
        return True
    except (ValueError, TypeError):
        return False
    

def parse_gc_percent(value):
    """
    Some fields in the gc_percent column of the refseq assembly file contain values over 100%.
    Helper function that gets rid of values over 100%.
    """
    try:
        val = float(value)
        if val < 0:  # extremely unlikely scenario but handling anyways
            return None
        elif val > 100:  # Cleanest way to handle values over 100
            return None
        else:
            return val
    except (ValueError, TypeError):
        return None


def iqr_filter(values):
    """
    Apply IQR filtering to a list of numeric values.
    """
    if len(values) < 4:
        return values

    q1 = np.percentile(values, 25)
    q3 = np.percentile(values, 75)
    iqr = q3 - q1

    lower = q1 - 1.5 * iqr
    upper = q3 + 1.5 * iqr

    return [v for v in values if lower <= v <= upper]


def summarize(values):
    if not values:
        return ["NA"] * 6

    return [
        min(values),
        max(values),
        statistics.median(values),
        statistics.mean(values),
        statistics.stdev(values) if len(values) >= 2 else 0,
        len(values)
    ]

def bp_to_mb(value):
    """
    Convert base pairs to megabase pairs (Mb).
    """
    return value / 1_000_000


def calculate_assembly_stats(assembly_summary_file):
    """
    Stream the RefSeq assembly summary file and compute
    per-taxid aggregated statistics.
    Output the NCBI_Assembly_stats_{YYMMDD}.txt file.
    """

    # Initialize the tax ID dictionary
    taxid_data = defaultdict(lambda: {
        "species": None,
        "genome_size": [],
        "gc_percent": [],
        "cds_count": []
    })

    logging.info("Reading assembly summary file...")


    with open_gzip(assembly_summary_file) as fh:
        for line in fh:
            if line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")

            try:
                taxid = int(fields[5])
                organism_name = fields[7]
                genome_size = fields[25]
                gc_percent = parse_gc_percent(fields[27])
                cds = fields[35]
            except (IndexError, ValueError):
                continue

            entry = taxid_data[taxid]

            if entry["species"] is None:
                entry["species"] = organism_name
            elif entry["species"] != organism_name:
                logging.debug(
                    f"TaxID {taxid} has multiple organism names: "
                    f"{entry['species']} vs {organism_name}"
                )

            if is_float(genome_size):
                entry["genome_size"].append(float(genome_size))

            if gc_percent is not None:
                entry["gc_percent"].append(gc_percent)

            if is_float(cds):
                entry["cds_count"].append(int(float(cds)))

    logging.info(f"Collected data for {len(taxid_data)} tax IDs")

    records = []

    for taxid, data in taxid_data.items():
        gs = iqr_filter(data["genome_size"])
        gc = iqr_filter(data["gc_percent"])
        cds = iqr_filter(data["cds_count"])

        gs_min, gs_max, gs_med, gs_mean, gs_sd, gs_n = summarize(gs)
        gc_min, gc_max, gc_med, gc_mean, gc_sd, gc_n = summarize(gc)
        cds_min, cds_max, cds_med, cds_mean, cds_sd, cds_n = summarize(cds)

        # Convert genome size stats to Mb
        if gs_min != "NA":
            gs_min = bp_to_mb(gs_min)
            gs_max = bp_to_mb(gs_max)
            gs_med = bp_to_mb(gs_med)
            gs_mean = bp_to_mb(gs_mean)
            gs_sd = bp_to_mb(gs_sd)

        records.append([
            data["species"],
            gs_min, gs_max, gs_med, gs_mean, gs_sd, gs_n,
            gc_min, gc_max, gc_med, gc_mean, gc_sd, gc_n,
            cds_min, cds_max, cds_med, cds_mean, cds_sd, cds_n,
            taxid
        ])

    return records


def main():
    args = parse_args()

    if args.version:
        print("assembly-stats-db v0.1.0")
        return
    
    records = calculate_assembly_stats(args.path_database)

    output_db = f"NCBI_Assembly_Stats_{timestamp}.txt"

    header = [
        "Species",
        "Min", "Max", "Median", "Mean", "StDev", "Assembly_count",
        "GC_Min", "GC_Max", "GC_Median", "GC_Mean", "GC_Stdev", "GC_count",
        "CDS_Min", "CDS_Max", "CDS_Median", "CDS_Mean", "CDS_Stdev", "CDS_count",
        "Consensus_TAXID"
    ]

    logging.info(f"Writing output file: {output_db}")

    with open(output_db, "w") as out:
        out.write("\t".join(header) + "\n")
        for record in records:
            out.write("\t".join(map(str, record)) + "\n")

if __name__ == "__main__":
    main()