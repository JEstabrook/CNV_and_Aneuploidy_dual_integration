
import os
import pandas as pd
import argparse

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
        df = pd.DataFrame([line.strip().strip().split('\t') for line in data],columns=header.split()[1:])  # Set column names excluding the '#h' prefix
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
        stitched_frame['diff'] = None
    else:
        cnv['diff'] = cnv.groupby(['chr', 'Type']).apply(calc_difference).T.values
        stitched_data = []
        cnv_calls = cnv.groupby(['chr', 'Type'])
        for (Chr, Type), cnv_frame in cnv_calls:
            stitched_cnv = collapse_rows(cnv_frame, cnv_stitch_window=cnv_stitch_window)
            stitched_data.append(stitched_cnv)
        stitched_frame = pd.concat(stitched_data)
    stitched_frame['Width'] = stitched_frame['End'] - stitched_frame['Start']
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

def calculate_overlap_percentage(set1, set2):
    """
    Calculates the overlap percentage of two sets of values.
    Arguments:
    set1 -- First set of values as a tuple (start1, end1).
    set2 -- Second set of values as a tuple (start2, end2).
    Returns:
    overlap_percentage -- The overlap percentage as a float.
    """
    start1, end1 = set1
    start2, end2 = set2
    # Calculate the overlap length
    overlap_length = min(end1, end2) - max(start1, start2)
    # If there's no overlap, return 0
    if overlap_length <= 0:
        return 0.0
    # Calculate the lengths of the two sets
    length1 = end1 - start1
    length2 = end2 - start2
    # Calculate the overlap percentage
    overlap_percentage = (overlap_length / min(length1, length2))
    return overlap_percentage

def pairwise_comparison(case, control, cnv_overlap_percentage=0.5, cnv_window=1000):
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
    control['unique_id'] = control[['Start','End','Width','diff']].sum(axis=1).map(hash)
    case['unique_id'] = case[['Start','End','Width','diff']].sum(axis=1).map(hash)
    grouped_case = case.groupby('chr')
    grouped_control = control.groupby('chr')
    filtered_rows = []
    case_overlap = []
    case_checked = []
    control_checked = []
    for chr_name, case_group in grouped_case:
        if chr_name in grouped_control.groups:
            control_group = grouped_control.get_group(chr_name)
            for _, case_row in case_group.iterrows():
                case_specific = case_row.copy().to_frame().T
                case_specific.columns = case_specific.columns + '_Case'
                case_id = case_specific['unique_id_Case'].values[0]
                n = 0
                for _, control_row in control_group.iterrows():
                    control_specific = control_row.copy().to_frame().T
                    control_specific.columns = control_specific.columns + '_Control'
                    control_id = control_specific['unique_id_Control'].values[0]
                    start_diff = abs(case_row['Start'] - control_row['Start'])
                    end_diff = abs(case_row['End'] - control_row['End'])
                    width_ratio = case_row['Width'] / control_row['Width']
                    overlap_percent = calculate_overlap_percentage(set([case_row['Start'],case_row['End']]),set([control_row['Start'],control_row['End']]))
                    if start_diff >= cnv_window and end_diff >= cnv_window and overlap_percent <= cnv_overlap_percentage and case_id not in case_checked and control_id not in control_checked:
                        filtered_rows.append(case_row)
                        break  # Move to the next case row
                    if n == 0:
                        if control_id not in control_checked:
                            case_control_overlap = pd.concat([case_specific.reset_index(drop=True),control_specific.reset_index(drop=True)],axis=1)
                            case_checked.append(case_id)
                            control_checked.append(control_id)
                            n +=1
                    else:
                        if (control_specific['End_Control'].values[0] > case_control_overlap['End_Control'].values[0]) and (control_specific['Start_Control'].values[0] < case_control_overlap['End_Case'].values[0]):
                            case_control_overlap['End_Control'] = control_specific['End_Control'].values[0]
                            control_checked.append(control_id)
                case_overlap.append(case_control_overlap)
        else:
            [filtered_rows.append(x[1].to_frame().T) for x in case_group.iterrows()]
    if len(filtered_rows) == 0:
        filtered_case = pd.DataFrame(filtered_rows,columns=case.columns)
    else:
        filtered_case = pd.concat(filtered_rows)
    if len(case_overlap) > 0:
        case_overlap_frame = pd.concat(case_overlap)
        case_overlap_frame['Width_Control'] = case_overlap_frame['End_Control'] - case_overlap_frame['Start_Control']
    else:
        case_overlap_frame = pd.DataFrame(columns=['Id_Case', 'chr_Case', 'Start_Case', 'End_Case', 'Width_Case', 'Type_Case', 'fractionalCopyNumber_Case', 'CopyNumber_Case', 'Confidence_Case', 'AlleleFreq_Case', 'Mask_overlap_fract_Case', 'Algorithm_Case', 'diff_Case', 'unique_id_Case', 'Id_Control', 'chr_Control', 'Start_Control', 'End_Control', 'Width_Control', 'Type_Control', 'fractionalCopyNumber_Control', 'CopyNumber_Control', 'Confidence_Control', 'AlleleFreq_Control', 'Mask_overlap_fract_Control', 'Algorithm_Control', 'diff_Control', 'unique_id_Control'])
    return filtered_case, case_overlap_frame

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
    merged = case.merge(control, on=['chr', 'types'], suffixes=('_case', '_control'), how='outer')
    filtered_rows = []
    merged_rows = []
    for _, row in merged.iterrows():
        if pd.isnull(row['fractChrLen_control']) or pd.isnull(row['fractCN_control']):
            filtered_rows.append(row)
        else:
            fractChrLen_diff = abs(row['fractChrLen_case'] - row['fractChrLen_control'])
            fractCN_diff = abs(row['fractCN_case'] - row['fractCN_control'])
            if fractChrLen_diff > aneuploidy_overlap_percentage and fractCN_diff > aneuploidy_overlap_percentage:
                filtered_rows.append(row)
            else:
                merged_rows.append(row)
    merged_case = pd.DataFrame(merged_rows,columns=merged.columns)
    filtered_case = pd.DataFrame(filtered_rows,columns=merged.columns).iloc[:,:5]
    filtered_case.columns = case.columns[:5]
    return filtered_case, merged_case

def process_cnvs(case, control, cnv_overlap_percentage=0.5, cnv_window=1000):
    """Processes CNV calls and reports unique and paired case & control CNV calls
    Args:
        case (pd.DataFrame): DataFrame representing the case data
        control (pd.DataFrame): DataFrame representing the control data
        cnv_overlap_percentage (float, optional): Maximum allowed reciprocal overlap percentage. Defaults to 0.5.
        cnv_window (int, optional): Window size for Start and End comparison. Defaults to 1000.
    Returns:
        pd.DataFrame: Filtered case DataFrame
    """
    cnv_reindex_cols = ['Cell type', 'Donor ID', 'Type', 'Id', 'Unique to case', 'chr', 'RefcontigID2', 'Start', 'End', 'Width', 'AlleleFreq', 'Case ID', 'Found in case', 'Self_molecule_count', 'Control ID', 'Found in control', 'Control_molecule_count', 'Algorithm', 'fractionalCopyNumber','Confidence']
    unique_case_cnvs, cnv_comp_table = pairwise_comparison(case=case,control=control, cnv_overlap_percentage=cnv_overlap_percentage, cnv_window=cnv_window)
    cnv_reindex_comp_cols = ['Cell type', 'Donor ID', 'Type_Case', 'Id_Case', 'Unique to case', 'chr_Case', 'RefcontigID2', 'Start_Case', 'End_Case', 'Width_Case', 'AlleleFreq_Case', 'Case ID', 'Found in case', 'Self_molecule_count', 'Control ID', 'Found in control', 'Control_molecule_count', 'Algorithm', 'fractionalCopyNumber','Start_Control', 'End_Control', 'Width_Control', 'AlleleFreq_Control','Confidence_Case','Confidence_Control']
    cnv_comp_table_indexed = cnv_comp_table.reindex(cnv_reindex_comp_cols, axis=1).rename(columns={'Type_Case':'Event Type', 'Id_Case':'Id','chr_Case':'ChrA location', 'RefcontigID2':'ChrB location','Start_Case':'Start','End_Case':'Stop', 'Width_Case':'SV size', 'Self_molecule_count':'Case count','Control_molecule_count':'Control count', 'AlleleFreq_Case':'VAF','AlleleFreq_Control':'VAF Control','Width_Control':'SV size Control','Confidence_Case':'Confidence','End_Control':'Stop_Control'})
    cnv_comp_table_indexed['CallType'] = ' CNV'
    cnv_comp_table_indexed['Algorithm'] = 'Region-based'
    cnv_comp_table_indexed['SV size control'] = cnv_comp_table_indexed['Stop_Control'] - cnv_comp_table_indexed['Start_Control']
    unique_case_cnvs_subset = unique_case_cnvs.reindex(cnv_reindex_cols, axis=1).rename(columns={'Type':'Event Type', 'chr':'ChrA location', 'RefcontigID2':'ChrB location','RefStartPos':'Start','End':'Stop', 'Width':'SV size', 'Self_molecule_count':'Case count','Control_molecule_count':'Control count', 'AlleleFreq':'VAF'})
    unique_case_cnvs_subset['CallType'] = ' CNV'
    return unique_case_cnvs_subset, cnv_comp_table_indexed

def process_aneuploidies(case, control, aneuploidy_overlap_percentage):
    """Processes Aneuploidy calls and reports unique and paired case & control Aneuploidy calls
    
    Args:
        case (pd.DataFrame): DataFrame representing the case data
        control (pd.DataFrame): DataFrame representing the control data
        aneuploidy_overlap_percentage (float, optional): Maximum allowed overlap percentage difference. Defaults to 0.5.
    
    Returns:
        unique_case_aneu_subset (pd.DataFrame): Rows unique to case DataFrame
        aneu_comp_table_indexed (pd.DataFrame): Aneuploidy events shared across control and case
    """
    contig_length = {1: 248956422.0, 2: 242193529.0, 3: 198295559.0, 4: 190214555.0, 5: 181538259.0, 6: 170805979.0, 7: 159345973.0, 8: 145138636.0, 9: 138394717.0, 10: 133797422.0, 11: 135086622.0, 12: 133275309.0, 13: 114364328.0, 14: 107043718.0, 15: 101991189.0, 16: 90338345.0, 17: 83257441.0, 18: 80373285.0, 19: 58617616.0, 20: 64444167.0, 21: 46709983.0, 22: 50818468.0, 23: 156040895.0, 24: 57227415.0}
    aneu_reindex_cols = ['Cell type', 'Donor ID', 'types', 'Id', 'Unique to case', 'chr', 'RefcontigID2', 'Start', 'End', 'fractChrLen', 'AlleleFreq', 'Case ID', 'Found in case', 'Self_molecule_count', 'Control ID', 'Found in control', 'Control_molecule_count', 'Algorithm', 'fractCN','score']
    if (case.shape[0] != 0) & (control.shape[0] !=0):
        unique_case_aneuploidies, aneu_comp_table = compare_aneuploidy_data(case=case, control=control, aneuploidy_overlap_percentage=aneuploidy_overlap_percentage)
        aneu_reindex_comp_cols = ['Cell type', 'Donor ID', 'types', 'Id', 'Unique to case', 'chr', 'RefcontigID2', 'Start', 'End', 'fractChrLen_case', 'AlleleFreq', 'Case ID', 'Found in case', 'Self_molecule_count', 'Control ID', 'Found in control', 'Control_molecule_count', 'Algorithm', 'fractCN_case', 'fractChrLen_control', 'fractCN_control', 'score_control','score_case']
        aneu_comp_table_indexed = aneu_comp_table.reindex(aneu_reindex_comp_cols, axis=1).rename(columns={'types':'Event Type', 'chr':'ChrA location', 'RefcontigID2':'ChrB location','RefStartPos':'Start','End':'Stop','Self_molecule_count':'Case count','Control_molecule_count':'Control count', 'AlleleFreq':'VAF','fractCN_case':'fractionalCopyNumber','fractCN_control':'fractionalCopyNumber_control','score_case':'Confidence', 'score_control':'Confidence_Control'})
        aneu_comp_table_indexed['CallType'] = 'Aneuploidy'
        aneu_comp_table_indexed['Algorithm'] = 'Region-based'
        aneu_comp_table_indexed['fractChrLen_case'] = aneu_comp_table_indexed.apply(lambda row: row['fractChrLen_case'] * contig_length[int(row['ChrA location'])], axis=1)
        aneu_comp_table_indexed['fractChrLen_control'] = aneu_comp_table_indexed.apply(lambda row: row['fractChrLen_control'] * contig_length[int(row['ChrA location'])], axis=1)
        aneu_comp_table_indexed = aneu_comp_table_indexed.rename(columns={'fractChrLen_case':'SV size','fractChrLen_control':'SV size Control'})
        unique_case_aneu_subset = unique_case_aneuploidies.reindex(aneu_reindex_cols, axis=1).rename(columns={'types':'Event Type', 'chr':'ChrA location', 'RefcontigID2':'ChrB location','RefStartPos':'Start','End':'Stop', 'fractChrLen':'SV size', 'Self_molecule_count':'Case count','Control_molecule_count':'Control count', 'AlleleFreq':'VAF','fractCN':'fractionalCopyNumber','score':'Confidence'})
        unique_case_aneu_subset['Algorithm'] = 'Region-based'
        unique_case_aneu_subset['CallType'] = 'Aneuploidy'
    else:
        empty_cols = ['Cell type', 'Donor ID', 'Event Type', 'Id', 'Unique to case', 'ChrA location', 'ChrB location', 'Start', 'Stop', 'SV size', 'VAF', 'Case ID', 'Found in case', 'Case count', 'Control ID', 'Found in control', 'Control count', 'Algorithm', 'fractionalCopyNumber', 'SV size control', 'fractionalCopyNumber_control', 'Confidence_control', 'Confidence', 'CallType']
        aneu_comp_table_indexed = pd.DataFrame(columns=empty_cols)
        unique_case_aneu_subset = pd.DataFrame(columns=aneu_reindex_cols)
    return unique_case_aneu_subset, aneu_comp_table_indexed

def process_smap(dual_smap_frame):
    """Processes SV calls and reports unique and paired case & control SV calls

    Args:
        dual_smap_frame (pd.DataFrame): DataFrame consisting of dual annotation SV calls
    """
    smap_reindex_cols = ['Cell type', 'Donor ID', 'Type', 'SmapEntryID', 'Unique to case', 'RefcontigID1', 'RefcontigID2', 'RefStartPos', 'RefEndPos', 'SVsize', 'VAF', 'Case ID', 'Found in case', 'Self_molecule_count', 'Control ID', 'Control_molecule_count', 'Algorithm', 'fractionalCopyNumber', 'Confidence','Found_in_control_sample_molecules','Found_in_control_sample_assembly']
    smap_subset = dual_smap_frame.reindex(smap_reindex_cols, axis=1).rename(columns={'SmapEntryID':'Id', 'Type':'Event Type', 'RefcontigID1':'ChrA location', 'RefcontigID2':'ChrB location','RefStartPos':'Start','RefEndPos':'Stop', 'SVsize':'SV size', 'Self_molecule_count':'Case count','Control_molecule_count':'Control count','Found_in_control_sample_molecules':'Found in control'})
    smap_subset['CallType'] = 'SV'
    smap_filtered = smap_subset[(smap_subset['Found in control'] != 'yes') & (smap_subset['Found_in_control_sample_assembly'] != 'yes')]
    return smap_filtered, smap_subset

def process_all_calls(smap_subset, aneu_comp_table_indexed, cnv_comp_table_indexed):
    """_summary_

    Args:
        smap_subset (_type_): _description_
        aneu_comp_table_indexed (_type_): _description_
        cnv_comp_table_indexed (_type_): _description_
    """
    reindex_cols = ['Cell type', 'Donor ID', 'Event Type', 'Id', 'ChrA location', 'ChrB location', 'Start','Start_Control', 'Stop', 'Stop_Control','SV size','SV size control','VAF', 'VAF Control', 'Case count', 'Control count','fractionalCopyNumber', 'fractionalCopyNumber_control', 'Confidence', 'Confidence_Control', 'Case ID', 'Control ID','Algorithm', 'Unique to case', 'Found in case', 'Found in control', 'Found_in_control_sample_assembly', 'CallType']
    smap_subset = smap_subset.reindex(reindex_cols,axis=1)
    aneu_comp_table_indexed = aneu_comp_table_indexed.reindex(reindex_cols,axis=1)
    cnv_comp_table_indexed = cnv_comp_table_indexed.reindex(reindex_cols,axis=1)
    smap_subset.reset_index(drop=True,inplace=True)
    aneu_comp_table_indexed.reset_index(drop=True,inplace=True)
    cnv_comp_table_indexed.reset_index(drop=True,inplace=True)
    all_calls = pd.concat([smap_subset, aneu_comp_table_indexed, cnv_comp_table_indexed],ignore_index=True,join='outer')
    convert_dict = {'VAF': float, 'VAF Control': float, 'Confidence':float, 'Confidence_Control':float}
    all_calls = all_calls.astype(convert_dict)
    return all_calls

def compare_calls(dual_aneuploidy, dual_smap, dual_cnv, control_aneuploidy, control_cnv, out_file, donor_id, case_id, celltype, cnv_overlap_percentage=0.3, aneuploidy_overlap_percentage=0.5, cnv_window=1000, cnv_stitch_window=550000):
    """This function compares case and control Aneuploidy and CNV calls and reports case specific SV, CNV and Aneuploidy calls

    Args:
        dual_aneuploidy (str): relative path to aneuploidy file from dual annotation zip
        dual_smap (str): relative path to smap file from dual annotation zip
        dual_cnv (str): relative path to CNV file from dual annotation zip
        control_aneuploidy (str): relative path to aneuploidy file from control zip
        control_cnv (str): relative path to CNV file from control zip
        out_file (str): outfile handle
        donor_id (str): donor id
        case_id (str): case id
        celltype (str): cell type
        cnv_overlap_percentage (float): maximum reciprocal overlap percentage to consider CNV unique to case
        aneuploidy_overlap_percentag (float): maximum fractional coverage difference to consider Aneuploidies unique to case
        cnv_window (int): base pair window to extend start and stop positions by for CNV calls
        cnv_stitch_window (int): base pair window used to extend and join neighboring cnvs that fall within n-basepairs from one another. Defaults to 550000.
    """
    dual_smap_frame = read_smap(dual_smap)
    dual_aneuploidy_frame = read_aneuploidy(dual_aneuploidy)
    dual_cnv_frame = read_cnv(dual_cnv)
    stitched_dual_cnvs = stitching(dual_cnv_frame, cnv_stitch_window=cnv_stitch_window)

    control_aneuploidy_frame = read_aneuploidy(control_aneuploidy)
    control_cnv_frame = read_cnv(control_cnv)
    stitched_control_cnvs = stitching(control_cnv_frame, cnv_stitch_window=cnv_stitch_window)

    unique_case_cnvs_subset, cnv_comp_table_indexed = process_cnvs(case=stitched_dual_cnvs, control=stitched_control_cnvs, cnv_overlap_percentage=cnv_overlap_percentage, cnv_window=cnv_window)
    unique_case_aneu_subset, aneu_comp_table_indexed = process_aneuploidies(case=dual_aneuploidy_frame, control=control_aneuploidy_frame, aneuploidy_overlap_percentage=aneuploidy_overlap_percentage)
    print(unique_case_aneu_subset)
    print(aneu_comp_table_indexed)

    smap_filtered, smap_subset = process_smap(dual_smap_frame)

    out_table = pd.concat([smap_filtered, unique_case_aneu_subset, unique_case_cnvs_subset],ignore_index=True)
    out_table['Case ID'] = case_id
    out_table['Cell type'] = celltype

    all_calls = process_all_calls(smap_subset, aneu_comp_table_indexed, cnv_comp_table_indexed)
    all_calls['Case ID'] = case_id
    all_calls['Cell type'] = celltype

    with pd.ExcelWriter(out_file) as writer:
            all_calls.to_excel(writer, sheet_name='DualAnnotationResults',index=False,float_format="%.3f")
            out_table.to_excel(writer, sheet_name='Case_Results',index=False,float_format="%.3f")

def main():
    parser = argparse.ArgumentParser(
        """This script compares Aneuploidy and CNV calls for samples with corresponding dual variant annotation results. This script will generate a table that contains case-specific SVs, CNVs and aneuploidy calls.""")
    parser.add_argument('--dual_aneuploidy', type=str, help="relative path to dual annotation *_Aneuploidy.txt")
    parser.add_argument('--dual_smap', type=str, help="relative path to dual annotation *_Annotated_SV.smap")
    parser.add_argument('--dual_cnv', type=str, help="relative path to dual annotation *_CNV.txt")
    parser.add_argument('--control_aneuploidy', type=str, help="relative path to control *_Aneuploidy.txt")
    parser.add_argument('--control_cnv', type=str, help="relative path to control *_CNV.txt")
    parser.add_argument('--out_file', type=str, help="output file handle")
    parser.add_argument('--case_id', type=str, help="Case ID")
    parser.add_argument('--celltype', type=str, help="cell type")
    parser.add_argument('--cnv_overlap_percentage', type=float, nargs='?', const=1, default=0.3, help="maximum reciprocal overlap ratio allowed for CNV calls to be considered unique")
    parser.add_argument('--aneuploidy_overlap_percentage', type=float, nargs='?', const=1, default=0.5,  help="maximum fractional difference allowed for Aneuploidy calls to be considered unique")
    parser.add_argument('--cnv_window', type=int, nargs='?', const=1, default=1000, help="bp window buffer at start and end of CNV calls")
    parser.add_argument('--cnv_stitch_window', type=int, nargs='?', const=1, default=550000, help="bp window to merge neighboring CNV calls")



    # parse command line arguments
    args = parser.parse_args()
    print(args)
    dual_aneuploidy = args.dual_aneuploidy
    dual_smap = args.dual_smap
    dual_cnv = args.dual_cnv
    control_aneuploidy = args.control_aneuploidy
    control_cnv = args.control_cnv
    cnv_window = args.cnv_window
    aneuploidy_overlap_percentage = args.aneuploidy_overlap_percentage
    cnv_overlap_percentage = args.cnv_overlap_percentage
    out_file = args.out_file
    donor_id = args.donor_id
    case_id = args.case_id
    celltype = args.celltype
    cnv_stitch_window = args.cnv_stitch_window
    
    compare_calls(dual_aneuploidy, dual_smap, dual_cnv, control_aneuploidy, control_cnv, out_file, donor_id, case_id, celltype, cnv_overlap_percentage, aneuploidy_overlap_percentage, cnv_window, cnv_stitch_window)


if __name__ == "__main__":
    main()