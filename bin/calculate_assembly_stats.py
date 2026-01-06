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


logging.basicConfig(level = logging.INFO, format = '%(levelname)s : %(message)s')

timestamp = datetime.now().strftime("%Y%m%d")

def parse_args(args=None):
    Description='Compare local assembly to expected assembly size based on taxonomy.'

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


def calculate_assembly_stats(url, target_taxid, sample_name, assembly_length, total_tax, sample_gc_percent, found):
    """
    Stream through the NCBI assembly_summary_refseq.txt file,
    compute mean genome_size (after IQR filtering) and mean gc_percent
    for all records matching the target_taxid.
    """

    logging.info(f"Fetching NCBI assembly summary for taxid {target_taxid} ...")

    genome_sizes = []
    gc_percents = []

    try:
        with urllib.request.urlopen(url) as response:
            for raw_line in response:
                line = raw_line.decode("utf-8").strip()
                if not line or line.startswith("#"):
                    continue

                cols = line.split("\t")
                if len(cols) <= 27:  # At minimum, column count should go to column 27 (zero-based numbering), which is 'gc_percent'. 
                    continue

                taxid = cols[5].strip()
                if taxid != str(target_taxid):
                    continue  # skip rows that are not target tax ID

                if total_tax == None:
                    total_tax = cols[7].strip()

                found = True

                genome_size_val = cols[25].strip()
                gc_percent_val = parse_gc_percent(cols[27].strip())  # Remove all gc_percent values over 100

                # Handle strings in genome_size and validate numeric fields
                if is_float(genome_size_val) and gc_percent_val is not None:
                    genome_size = float(genome_size_val)
                    gc_percent = float(gc_percent_val)
                    genome_sizes.append(genome_size)
                    gc_percents.append(gc_percent)
                else:
                    logging.debug(
                        f"Skipping malformed entry for taxid {taxid}: genome_size='{genome_size_val}', gc_percent='{gc_percent_val}'"
                    )
        
        if not found:
            logging.warning(f"No NCBI matches found for target taxid '{target_taxid}'")
            return None

        if not genome_sizes:
            logging.warning(f"No valid genome entries found for taxid {target_taxid}")
            return None, None

        # --- IQR filtering on genome_size and gc_percent ---
        Q1 = np.percentile(genome_sizes, 25)
        Q3 = np.percentile(genome_sizes, 75)
        IQR = Q3 - Q1

        lower_bound = Q1 - 1.5 * IQR
        upper_bound = Q3 + 1.5 * IQR

        filtered_sizes = [x for x in genome_sizes if lower_bound <= x <= upper_bound]  # Remove outliers

        Q1_gc = np.percentile(gc_percents, 25)
        Q3_gc = np.percentile(gc_percents, 75)
        IQR_gc = Q3_gc - Q1_gc

        gc_lower_bound = Q1_gc - 1.5 * IQR_gc
        gc_upper_bound = Q3_gc + 1.5 * IQR_gc

        filtered_gc = [x for x in gc_percents if gc_lower_bound <= x <= gc_upper_bound]  # Remove GC% outliers

        if not filtered_sizes or not filtered_gc:
            logging.warning(f"All values filtered out for taxid {target_taxid}")
            return None

        # Final stats
        expected_length = statistics.mean(filtered_sizes)  # Mean genome size
        # For continuity, only calculate std dev for species with 10 or more references
        if len(filtered_sizes) >= 10:
            stdev_genome_size = statistics.stdev(filtered_sizes)
            if int(assembly_length) > int(expected_length):
                bigger = int(assembly_length)
                smaller = int(expected_length)
            else:
                smaller = int(assembly_length)
                bigger = int(expected_length)

            logging.debug("Calculating the standard deviations")
            stdevs = (bigger - smaller) / stdev_genome_size

        else:
            stdev_genome_size = "Not calculated on species with n<10 references"
            stdevs = "NA"
        
        ## GC
        species_gc_mean = statistics.mean(gc_percents)
        gc_min = min(gc_percents)
        gc_max = max(gc_percents)
        gc_count = len(gc_percents)
        species_gc_percent_stdev = statistics.stdev(filtered_gc) if len(filtered_gc) >= 10 else "Not calculated on species with n<10 references" 


        logging.info(
            f"Taxid {target_taxid}: mean genome_size (IQR-filtered) = {expected_length:.2f}, "
            f"mean GC% = {species_gc_mean:.2f} (n={len(filtered_sizes)})"
        )

        # with open(f"NCBI_Assembly_stats_{timestamp}.txt", 'w') as outfile:
        #             outfile.write(f"Sample: {sample_name}\nTax: {total_tax}\nNCBI_TAXID: {target_taxid}\nSpecies_GC_StDev: {species_gc_percent_stdev}\nSpecies_GC_Min: {gc_min}\nSpecies_GC_Max: {gc_max}\nSpecies_GC_Mean: {species_gc_mean}\nSpecies_GC_Count: {gc_count}\nSample_GC_Percent: {sample_gc_percent}")

        # List of column headers: Species, Min, Max, Median, Mean, StDev, Assembly_count, GC_Min, GC_Max, GC_Median, GC_Mean, GC_Stdev, GC_count, CDS_Min, CDS_Max, CDS_Median, CDS_Mean, CDS_Stdev, CDS_count, Consensus_TAXID
    
    except Exception as e:
        logging.error(f"Error computing taxid genome stats for {target_taxid}: {e}")
        return None