#!/usr/bin/env python3
import argparse
import re
import os
import logging
import sys
import pandas as pd
import urllib.request
import numpy as np
import statistics
from datetime import datetime
import gzip


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


def calculate_assembly_stats(assembly_summary_file):
    """
    Process the NCBI assembly_summary_refseq.txt file,
    compute mean genome_size (after IQR filtering) and mean gc_percent
    for all records matching the target_taxid.
    Output the NCBI_Assembly_stats_{YYMMDD}.txt file.
    """


    try:
        with gzip.open(assembly_summary_file, "rt") as asf:
                for line in asf:
                    if not line or line.startswith("#"):
                        continue


        # with open(f"NCBI_Assembly_stats_{timestamp}.txt", 'w') as outfile:
        #             outfile.write(f"Sample: {sample_name}\nTax: {total_tax}\nNCBI_TAXID: {target_taxid}\nSpecies_GC_StDev: {species_gc_percent_stdev}\nSpecies_GC_Min: {gc_min}\nSpecies_GC_Max: {gc_max}\nSpecies_GC_Mean: {species_gc_mean}\nSpecies_GC_Count: {gc_count}\nSample_GC_Percent: {sample_gc_percent}")

        # List of column headers: Species, Min, Max, Median, Mean, StDev, Assembly_count, GC_Min, GC_Max, GC_Median, GC_Mean, GC_Stdev, GC_count, CDS_Min, CDS_Max, CDS_Median, CDS_Mean, CDS_Stdev, CDS_count, Consensus_TAXID
    
    except Exception as e:
        logging.error(f"Error computing taxid genome stats for {target_taxid}: {e}")
        return None