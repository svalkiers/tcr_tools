import pandas as pd
import os

from parsing_functions import select_highest_duplicate_count

def read_10x(path):
    """
    Read and parse 10X Cellranger output. Note that this function removes any
    cells for which only a single chain was detected. Double alpha and or beta
    chains are also removed.
    """
    cols = ['cell_id', 'clone_id', 'sequence_id', 'sequence', 'productive', 'v_call', 'j_call', 'junction', 'junction_aa', 'duplicate_count']
    df = pd.read_csv(path, sep="\t")

    df = df[cols]
    df = df[df.productive=="T"]
    df = df[df.v_call.str.contains(("TRA|TRB"))] # keep alphabeta only
    df["chain"] = df.v_call.apply(lambda x: "B" if "TRBV" in x else "A")

    # Select one alpha-beta pair per cell based on duplicate count
    df = df.groupby('cell_id').apply(select_highest_duplicate_count).reset_index(drop=True)

    return df

def demultiplex_10x_paired_tcr(vdj_file, barcode_file, write_to_file=False):
    """
    Split 10X VDJ file into subfiles per sample. This function requires a barcode file that
    must be manually created. The barcode file contains two columns, one carrying information
    about the cell barcodes and the other the different sample tags.
    """
    df = read_10x(vdj_file)
    barcodes = pd.read_csv(barcode_file)
    df = df.merge(barcodes, on="cell_id")
    clones = df.groupby("clone_id").sample_id.unique().index
    conditions = df.groupby("clone_id").sample_id.unique().values
    fp = []
    for idx, i in zip(clones, conditions):
        if len(i) > 1:
            cond = set([j[4:7] for j in i])
            if len(cond) == 2:
                if 'poo' not in cond:
                    fp.append(idx)
            elif len(cond) > 2:
                fp.append(idx)
    df = df[~df.clone_id.isin(fp)]

    if write_to_file:
        directory = os.path.dirname(barcode_file)
        out_folder = os.path.join(directory, "samples")
        if not os.path.exists(out_folder):
            os.mkdir(out_folder)
        samples = df.sample_id.unique()
        for sample in samples:
            print(sample)
            df[df["sample_id"]==sample].to_csv(f"{out_folder}/{sample}.tsv", sep="\t", index=False)
            
    return df

def to_paired_format(df):
    """
    Reformat DataFrame to pair alpha and beta chain information into a single row.
    """
    df = df.pivot(index="cell_id", columns="chain")
    df.columns = ['_'.join(i) for i in df.columns]

    cols_to_keep = {
        'clone_id_A' : 'clone_id', 
        'sequence_id_A' : 'sequence_a_id', 
        'sequence_id_B' : 'sequence_b_id', 
        'sequence_A' : 'sequence_a', 
        'sequence_B' : 'sequence_b', 
        'productive_A' : 'productive', 
        'v_call_A' : 'v_a_call', 
        'v_call_B' : 'v_b_call', 
        'j_call_A' : 'j_a_call', 
        'j_call_B' : 'j_b_call', 
        'junction_A' : 'junction_a', 
        'junction_B' : 'junction_b', 
        'junction_aa_A' : 'junction_a_aa', 
        'junction_aa_B' : 'junction_b_aa', 
        'duplicate_count_A' : 'duplicate_count_a',
        'duplicate_count_B' : 'duplicate_count_b'
        }

    if "sample_id_A" in df.columns:
        cols_to_keep["sample_id_A"] = "sample_id"
    else:
        pass

    df = df[list(cols_to_keep.keys())]
    df = df.rename(columns = cols_to_keep)

    return df