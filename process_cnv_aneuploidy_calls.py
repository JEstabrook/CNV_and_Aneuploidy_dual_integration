
import os
import pandas as pd
import argparse
import numpy as np
import shutil


def read_smap(file_path):
    """ Function parses smap file

    Args:
        file_path (str): relative path to smap file to be parsed

    Returns:
        df (pd.DataFrame) : DataFrame of SV calls extracted from smap file
    """
    with open(file_path, 'r') as file:
        lines = file.readlines()
    start_index = None
    for i, line in enumerate(lines):
        if line.startswith("#h"):
            start_index = i
            break
    if start_index is not None:
        header = lines[start_index].strip()  # Store the header line
        data = lines[start_index + 2:]  # Skip the header line and the following line
        df = pd.DataFrame([line.strip().split() for line in data],columns=header.split()[1:])  # Set column names excluding the '#h' prefix
        return df
    else:
        return None  # Return None if no header line is found


def read_aneuploidy(file_path):
    """Function parses aneuploidy calls

    Args:
        file_path (str): relative path to aneuploidy file to be parsed

    Returns:
        df (pd.DataFrame): DataFrame of Aneuploidy calls extracted from aneuploidy file
    """
    with open(file_path, 'r') as file:
        lines = file.readlines()
    start_index = None
    for i, line in enumerate(lines):
        if line.startswith("#chr"):
            start_index = i
            break
    if start_index is not None:
        header = lines[start_index].strip()  # Store the header line
        data = lines[start_index + 1:]  # Skip only the header line
        df = pd.DataFrame([line.strip().split() for line in data],columns=header.split()).rename(columns={'#chr':'chr'}) # Set column names
        return df
    else:
        return None 

def read_cnv(file_path):
    """Function parses CNV calls

    Args:
        file_path (str): relative path to CNV file to be parsed

    Returns:
        df (pd.DataFrame): DataFrame of CNV calls extracted from CNV file
    """
    with open(file_path, 'r') as file:
        lines = file.readlines()
    start_index = None
    for i, line in enumerate(lines):
        if line.startswith("#Id"):
            start_index = i
            break
    if start_index is not None:
        header = lines[start_index].strip()  # Store the header line
        data = lines[start_index + 1:]  # Skip only the header line
        df = pd.DataFrame([line.strip().split() for line in data],columns=header.split()).rename(columns={'#Id':'Id','Chromosome':'chr'}) # Set column names
        df['Start'] = df['Start'].astype(int)
        df['End'] = df['End'].astype(int)
        return df
    else:
        return None 

def calc_difference(cnv):
    """Calculate the difference between start and end positions of subsequent rows based on chromosome and SVtype

    Args:
        cnv (pd.DataFrame): DataFrame consisting of CNV calls

    Returns:
        pd.Series: base pair difference between CNV segments on the same chromosome
    """
    vals = (cnv['Start'] - cnv['End'].shift(1)).fillna(cnv.groupby(['chr', 'Type'])['Start'].transform('first'))
    return vals

def stitching(cnv, cnv_stitch_window=550000, width=True):
    """Stitchs CNV calls within a particular window bin together

    Args:
        cnv (pd.DataFrame): DataFrame of CNV calls extracted from CNV file parsed by read_cnv()
        cnv_stitch_window (int, optional): Designated window to collapse CNV calls by. Defaults to 550000.
        width (bool, optional): Return the width of the CNV segment. Defaults to True.

    Returns:
        pd.DataFrame: CNV stitched calls
    """
    cnv = cnv.sort_values(by=['chr', 'Type', 'Start'])
    if cnv.shape[0] == 0:
        stitched_frame = cnv
    else:
        cnv['diff'] = cnv.groupby(['chr', 'Type']).apply(calc_difference).T.values
        stitched_data = []
        cnv_calls = cnv.groupby(['chr', 'Type'])
        for (Chr, Type), cnv_frame in cnv_calls:
            stitched_cnv = collapse_rows(cnv_frame, cnv_stitch_window=cnv_stitch_window)
            stitched_data.append(stitched_cnv)
        stitched_frame = pd.concat(stitched_data)
    return stitched_frame

def collapse_rows(df, cnv_stitch_window):
    """Collapses rows in a DataFrame based on a condition
    Args:
        df (pd.DataFrame): Input DataFrame
        cnv_stitch_window (int): Threshold value for collapsing rows
    Returns:
        pd.DataFrame: Collapsed DataFrame
    """
    collapsed_data = []
    current_start = None
    current_end = None
    current_row = None
    for index, row in df.iterrows():
        if current_start is None:
            current_start = row['Start']
            current_end = row['End']
            current_row = row
        else:
            diff = row['diff']
            if diff < cnv_stitch_window:
                current_end = row['End']
                current_row['End'] = current_end  # Update the 'End' value in the current row
            else:
                collapsed_data.append(current_row)
                current_start = row['Start']
                current_end = row['End']
                current_row = row
    if current_row is not None:
        collapsed_data.append(current_row)
    collapsed_df = pd.DataFrame(collapsed_data)
    return collapsed_df

def pairwise_comparison(case, control, cnv_overlap_percentage, cnv_window):
    """Performs pairwise comparison between case and control DataFrames
    Args:
        case (pd.DataFrame): DataFrame representing the case data
        control (pd.DataFrame): DataFrame representing the control data
        cnv_overlap_percentage (float, optional): Maximum allowed reciprocal overlap percentage. Defaults to 0.5.
        cnv_window (int, optional): Window size for Start and End comparison. Defaults to 1000.
    Returns:
        pd.DataFrame: Filtered case DataFrame
    """
    case['Width'] = case['Width'].astype(float)
    control['Width'] = control['Width'].astype(float)
    filtered_rows = []
    for _, case_row in case.iterrows():
        for _, control_row in control.iterrows():
            start_diff = abs(case_row['Start'] - control_row['Start'])
            end_diff = abs(case_row['End'] - control_row['End'])
            width_ratio = case_row['Width'] / control_row['Width']
            if start_diff <= cnv_window and end_diff <= cnv_window and width_ratio <= cnv_overlap_percentage:
                filtered_rows.append(case_row)
                break  # Move to the next case row
    filtered_case = pd.DataFrame(filtered_rows, columns=case.columns)
    return filtered_case

def compare_aneuploidy_data(case, control, aneuploidy_overlap_percentage):
    """Compares case and control DataFrames and returns rows unique to case based on specified conditions
    
    Args:
        case (pd.DataFrame): DataFrame representing the case data
        control (pd.DataFrame): DataFrame representing the control data
        aneuploidy_overlap_percentage (float, optional): Maximum allowed overlap percentage difference. Defaults to 0.5.
    
    Returns:
        pd.DataFrame: Rows unique to case DataFrame
    """
    case['fractChrLen'] = case['fractChrLen'].astype(float)
    control['fractChrLen'] = control['fractChrLen'].astype(float)
    case['fractCN'] = case['fractCN'].astype(float)
    control['fractCN'] = control['fractCN'].astype(float)
    merged = case.merge(control, on=['chr', 'types'], suffixes=('_case', '_control'))
    filtered_rows = []
    for _, row in merged.iterrows():
        fractChrLen_diff = abs(row['fractChrLen_case'] - row['fractChrLen_control'])
        fractCN_diff = abs(row['fractCN_case'] - row['fractCN_control'])
        if fractChrLen_diff > aneuploidy_overlap_percentage and fractCN_diff > aneuploidy_overlap_percentage:
            filtered_rows.append(row)
    filtered_case = pd.DataFrame(filtered_rows, columns=case.columns)
    return filtered_case

def compare_calls(dual_aneuploidy, dual_smap, dual_cnv, control_aneuploidy, control_cnv, out_file, cnv_overlap_percentage=0.3, aneuploidy_overlap_percentage=0.5, cnv_window=1000, cnv_stitch_window=550000):
    """This function compares case and control Aneuploidy and CNV calls and reports case specific SV, CNV and Aneuploidy calls

    Args:
        dual_aneuploidy (str): relative path to aneuploidy file from dual annotation zip
        dual_smap (str): relative path to smap file from dual annotation zip
        dual_cnv (str): relative path to CNV file from dual annotation zip
        control_aneuploidy (str): relative path to aneuploidy file from control zip
        control_cnv (str): relative path to CNV file from control zip
        out_file (str): outfile handle
        cnv_overlap_percentage (float): maximum reciprocal overlap percentage to consider CNV unique to case
        aneuploidy_overlap_percentag (float): maximum fractional coverage difference to consider Aneuploidies unique to case
    """

    dual_smap_frame = read_smap(dual_smap)
    dual_aneuploidy_frame = read_aneuploidy(dual_aneuploidy)
    dual_cnv_frame = read_cnv(dual_cnv)
    stitched_dual_cnvs = stitching(dual_cnv_frame, cnv_stitch_window=cnv_stitch_window)

    control_aneuploidy_frame = read_aneuploidy(control_aneuploidy)
    control_cnv_frame = read_cnv(control_cnv)
    stitched_control_cnvs = stitching(control_cnv_frame, cnv_stitch_window=cnv_stitch_window)

    cnv_reindex_cols = ['Cell type', 'Donor ID', 'Type', 'Id', 'Unique to case', 'chr', 'RefcontigID2', 'Start', 'End', 'Width', 'AlleleFreq', 'Case ID', 'Found in case', 'Self_molecule_count', 'Control ID', 'Found in control', 'Control_molecule_count', 'Algorithm', 'fractionalCopyNumber','Confidence']
    unique_case_cnvs = pairwise_comparison(case=stitched_dual_cnvs,control=stitched_control_cnvs, cnv_overlap_percentage=cnv_overlap_percentage, cnv_window=cnv_window)
    unique_case_cnvs_subset = unique_case_cnvs.reindex(cnv_reindex_cols, axis=1).rename(columns={'Type':'Event Type', 'chr':'Chr1 location', 'RefcontigID2':'Chr2 location','RefStartPos':'Start','RefEndPos':'Stop', 'Width':'SV size', 'Self_molecule_count':'Case count','Control_molecule_count':'Control count', 'AlleleFreq':'VAF'})
    unique_case_cnvs_subset['CallType'] = ' CNV'

    aneu_reindex_cols = ['Cell type', 'Donor ID', 'types', 'Id', 'Unique to case', 'chr', 'RefcontigID2', 'Start', 'End', 'fractChrLen', 'AlleleFreq', 'Case ID', 'Found in case', 'Self_molecule_count', 'Control ID', 'Found in control', 'Control_molecule_count', 'Algorithm', 'fractCN','score']
    unique_case_aneuploidies = compare_aneuploidy_data(case=dual_aneuploidy_frame, control=control_aneuploidy_frame, aneuploidy_overlap_percentage=aneuploidy_overlap_percentage)
    unique_case_aneu_subset = unique_case_aneuploidies.reindex(aneu_reindex_cols, axis=1).rename(columns={'types':'Event Type', 'chr':'Chr1 location', 'RefcontigID2':'Chr2 location','RefStartPos':'Start','RefEndPos':'Stop', 'fractChrLen':'SV size', 'Self_molecule_count':'Case count','Control_molecule_count':'Control count', 'AlleleFreq':'VAF','fractCN':'fractionalCopyNumber'})
    unique_case_aneu_subset['CallType'] = 'Aneuploidy'

    smap_reindex_cols = ['Cell type', 'Donor ID', 'SVType', 'SmapEntryID', 'Unique to case', 'RefcontigID1', 'RefcontigID2', 'RefStartPos', 'RefEndPos', 'SVsize', 'VAF', 'Case ID', 'Found in case', 'Self_molecule_count', 'Control ID', 'Found in control', 'Control_molecule_count', 'Algorithm', 'fractionalCopyNumber', 'Confidence']
    smap_subset = dual_smap_frame.reindex(smap_reindex_cols, axis=1).rename(columns={'SmapEntryID':'Donor ID', 'SVType':'Event Type', 'RefcontigID1':'Chr1 location', 'RefcontigID2':'Chr2 location','RefStartPos':'Start','RefEndPos':'Stop', 'SVsize':'SV size', 'Self_molecule_count':'Case count','Control_molecule_count':'Control count'})
    smap_subset['CallType'] = 'SV'

    out_table = pd.concat([smap_subset, unique_case_aneu_subset, unique_case_cnvs_subset],axis=1)
    with pd.ExcelWriter(out_file) as writer:
            out_table.to_excel(writer, sheet_name='DualAnnotationResults',index=False)

def main():
    parser = argparse.ArgumentParser(
        """This script compares Aneuploidy and CNV calls for samples with corresponding dual variant annotation results. This script will generate a table that contains case-specific SVs, CNVs and aneuploidy calls.""")
    parser.add_argument('--dual_aneuploidy', type=str, help="relative path to dual annotation *_Aneuploidy.txt")
    parser.add_argument('--dual_smap', type=str, help="relative path to dual annotation *_Annotated_SV.smap")
    parser.add_argument('--dual_cnv', type=str, help="relative path to dual annotation *_CNV.txt")
    parser.add_argument('--control_aneuploidy', type=str, help="relative path to control *_Aneuploidy.txt")
    parser.add_argument('--control_cnv', type=str, help="relative path to control *_CNV.txt")
    parser.add_argument('--out_file', type=str, help="output file handle")
    parser.add_argument('--cnv_overlap_percentage', type=float, nargs='?', const=1, default=0.3, help="maximum reciprocal overlap ratio allowed for CNV calls to be considered unique")
    parser.add_argument('--aneuploidy_overlap_percentage', type=float, nargs='?', const=1, default=0.5,  help="maximum fractional difference allowed for Aneuploidy calls to be considered unique")
    parser.add_argument('--cnv_window', type=int, nargs='?', const=1, default=1000, help="bp window buffer at start and end of CNV calls")
    parser.add_argument('--cnv_stitch_window', type=int, nargs='?', const=1, default=550000, help="bp window to merge neighboring CNV calls")



    # parse command line arguments
    args = parser.parse_args()
    dual_aneuploidy = args.dual_aneuploidy
    dual_smap = args.dual_smap
    dual_cnv = args.dual_cnv
    control_aneuploidy = args.control_aneuploidy
    control_cnv = args.control_cnv
    cnv_window = args.cnv_window
    aneuploidy_overlap_percentage = args.aneuploidy_overlap_percentage
    cnv_overlap_percentage = args.cnv_overlap_percentage
    out_file = args.out_file
    cnv_stitch_window = args.cnv_stitch_window
    
    compare_calls(dual_aneuploidy, dual_smap, dual_cnv, control_aneuploidy, control_cnv, out_file, cnv_overlap_percentage, aneuploidy_overlap_percentage, cnv_window, cnv_stitch_window)


if __name__ == "__main__":
    main()