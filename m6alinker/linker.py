#!/usr/bin/env python

import pandas as pd
import numpy as np
import gffutils
from collections import defaultdict
import warnings

def process_files(input_file, gtf_file, output_prefix):
    print(f"Loading input site file: {input_file}")
    print(f"Loading GTF file: {gtf_file}")

    # Load modification sites
    modification_sites = pd.read_csv(input_file)

    # Create or load GTF database
    db_file = gtf_file + ".db"
    try: # Check whether the database already exists
        db = gffutils.FeatureDB(db_file, keep_order=True)
    except:
        print("Creating GTF database...")
        db = gffutils.create_db(gtf_file, dbfn=db_file, force=True, keep_order=True, disable_infer_genes=True, disable_infer_transcripts=True)
        db = gffutils.FeatureDB(db_file, keep_order=True)

    # Extract transcript IDs from the modification sites
    transcript_ids_to_keep = set(modification_sites["transcript_id"])
    
    # Filter features from the GTF database to include only relevant transcript IDs
    transcript_data = []
    filtered_features = []
    for feature in db.all_features():
        if "transcript_id" in feature.attributes:
            transcript_id = feature.attributes["transcript_id"][0].split(".")[0]
            if transcript_id in transcript_ids_to_keep:
                filtered_features.append(feature)  # Keep only relevant features
                if feature.featuretype == "transcript":
                    transcript_data.append({
                        "transcript_id": feature.id.split(".")[0],
                        "gene_name": feature.attributes.get("gene_name", [None])[0],
                        "transcript_type": feature.attributes.get("transcript_type", [None])[0],
                    })
    
    # Turn the transcrip data to DataFrame and drop dupilcated data
    transcript_data = pd.DataFrame(transcript_data).drop_duplicates("transcript_id")

    
    # Group filtered features by transcript_id
    grouped_features = defaultdict(list)
    for feature in filtered_features:
        transcript_id = feature.attributes["transcript_id"][0].split(".")[0]  # Remove version
        grouped_features[transcript_id].append(feature)

    # Get the transcripts features
    tx_features = calculate_transcript_features(grouped_features, transcript_data)
    
    # Filter exons for relevant transcripts
    exons = defaultdict(list)
    for tx_id, features in grouped_features.items():
        for feature in features:
            if feature.featuretype == "exon":
                exons[tx_id].append(feature)
    
    # Apply the convert_to_genome_coordinates function to each row of the modification_sites
    modification_sites["annotation"] = modification_sites.apply(
        lambda row: convert_to_genome_coordinates(row["transcript_id"], row["transcript_position"],exons, tx_features), axis=1
    )
    # Expand dictionary columns
    annotation_df = modification_sites["annotation"].apply(pd.Series)
    modification_sites = pd.concat([modification_sites, annotation_df], axis=1).drop(columns=["annotation"])

    # Merge with transcript biotype information
    annotated_modification_sites = modification_sites.merge(transcript_data, on="transcript_id", how="left")
    
    # Sort the annotated result according to the transcript ids
    annotated_modification_sites = annotated_modification_sites.sort_values(by='transcript_id')

    # Save output
    output_file = f"{output_prefix}_{input_file.split('/')[-1]}"
    annotated_modification_sites.to_csv(output_file, index=False)
    print(f"Annotated modification sites saved to {output_file}")

def calculate_transcript_features(grouped_features, transcript_data):
    '''
    Compute transcript features for grouped transcript data.
    For each transcript:
    - Sort CDS features by start and end positions to find start_codon and stop_codon positions.
    - Calculate CDS length by summing all CDS feature lengths.
    - Determine UTR5 and UTR3 lengths based on UTR positions relative to start_codon and stop_codon, considering the strand.
    - If no CDS features exist, set CDS_len, utr5_len, and utr3_len to 0.
    '''

    transcript_features = []  # Store the results for all transcripts

    for tx_id, features in grouped_features.items():
        cds_length = 0
        utr5_length = 0
        utr3_length = 0
        tx_length = 0
        start_codon_pos = None
        stop_codon_pos = None

        # Separate features by type
        cds_features = [f for f in features if f.featuretype == "CDS"]
        utr_features = [f for f in features if f.featuretype == "UTR"]
        exon_features = [f for f in features if f.featuretype == "exon"]


        # If no CDS features exist, set all lengths to 0， and compute transcript length by adding exons
        if not cds_features:
            cds_length = 0
            utr5_length = 0
            utr3_length = 0
            for exon in exon_features:
                tx_length += exon.end - exon.start + 1

        else:
            # Sort CDS features by start and end positions
            sorted_cds = sorted(cds_features, key=lambda x: (x.start, x.end))
            if utr_features[0].strand == "+":  # Positive strand
                start_codon_pos = sorted_cds[0].start  # First CDS start
                stop_codon_pos = sorted_cds[-1].end  # Last CDS end
            elif utr_features[0].strand == "-":  # Negative strand
                start_codon_pos = sorted_cds[-1].end  # Last CDS end
                stop_codon_pos = sorted_cds[0].start  # First CDS start

            # Calculate CDS length by summing the lengths of all CDS features
            for cds in sorted_cds:
                cds_length += cds.end - cds.start + 1

            # Handle UTRs and classify them as UTR5 or UTR3
            if utr_features:
                strand = utr_features[0].strand  # Assume all UTRs share the same strand
                for utr in utr_features:
                    utr_length = utr.end - utr.start + 1
                    if strand == "+":
                        if utr.end < start_codon_pos:  # 5' UTR is before start_codon
                            utr5_length += utr_length
                        elif utr.start > stop_codon_pos:  # 3' UTR is after stop_codon
                            utr3_length += utr_length
                    elif strand == "-":
                        if utr.start > start_codon_pos:  # 5' UTR is after start_codon on negative strand
                            utr5_length += utr_length
                        elif utr.end < stop_codon_pos:  # 3' UTR is before stop_codon on negative strand
                            utr3_length += utr_length
                
            tx_length = cds_length + utr3_length + utr5_length
    
        # Add the results for this transcript
        transcript_features.append({
            "transcript_id": tx_id,
            "tx_length": tx_length,
            "cds_length": cds_length,
            "utr5_length": utr5_length,
            "utr3_length": utr3_length
        })

    # Turn the transcript features into DataFrame and merge it with the 
    tx_features = pd.DataFrame(transcript_features)
    tx_features = tx_features.merge(transcript_data, on="transcript_id", how="left")

    return tx_features

# Function to convert transcript coordinates to genomic coordinates
def convert_to_genome_coordinates(tx_name, tx_pos, exonsdb, txfdb):
    """
    Convert transcript coordinates to genome coordinates and annotate regions.
    """
    if tx_name not in exonsdb:
        warnings.warn(f"Transcript {tx_name} not found in exons database.")
        return {
            "genome_pos": None,
            "chromosome": None,
            "dist_up_exon_junc": None,
            "dist_down_exon_junc": None,
            "region": None
        }
    
    tx_exons = sorted(exonsdb[tx_name], key=lambda x: x.start)  # Sort exons by start
    chrom = tx_exons[0].chrom
    strand = tx_exons[0].strand
    
    # Calculate cumulative exon lengths
    exon_lengths = [exon.end - exon.start + 1 for exon in tx_exons]
    cum_lengths = np.cumsum(exon_lengths)

    # Identify the exon containing the position
    idx = np.searchsorted(cum_lengths, tx_pos)
    
    if idx >= len(tx_exons):
        warnings.warn(f"Transcript position {tx_pos} is out of range for transcript {tx_name}")
        return {
            "genome_pos": None,
            "chromosome": None,
            "dist_up_exon_junc": None,
            "dist_down_exon_junc": None,
            "region": None
        }
    
    # Compute genomic position
    rel_pos = tx_pos - (cum_lengths[idx - 1] if idx > 0 else 0)
    exon = tx_exons[idx]
    if strand == "+":
        genome_pos = exon.start + rel_pos - 1
    else:
        genome_pos = exon.end - rel_pos + 1
    
    # Determine distance to exon junctions
    dist_upstream = rel_pos
    dist_downstream = exon.end - exon.start - rel_pos

    # Determine region (UTR, CDS, etc.)
    tx_feat = txfdb[txfdb["transcript_id"] == tx_name].iloc[0]
    if tx_feat["cds_length"] == 0:
        region = "ncRNA"
    elif tx_pos <= tx_feat["utr5_length"]:
        region = "UTR5"
    elif tx_pos <= tx_feat["utr5_length"] + tx_feat["cds_length"]:
        region = "CDS"
    elif tx_pos <= tx_feat["tx_length"]:
        region = "UTR3"
    else:
        region = None
    
    return {
        "genome_pos": genome_pos,
        "chromosome": chrom,
        "dist_up_exon_junc": dist_upstream,
        "dist_down_exon_junc": dist_downstream,
        "region": region
    }
