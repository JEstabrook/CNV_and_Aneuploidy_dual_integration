## Bionano - CNV and Aneuploidy dual annotation call integration
---
This script is designed to integrate and compare SV, CNV and Aneuploidy calls from case and control samples with paired dual variant annotation results. This script will generate an .xlsx that contains
SVs, CNVs and Aneuploidies unique to the designated case.

## SETUP
---
Current implementation of this script utilizes the following packages
```
Python (3.10.4)
Pandas (1.4.3)
argparse (1.1)
numpy (1.23.2)
docx (0.8.11)
```
## Using Conda

**Conda installation**

```
wget https://repo.anaconda.com/miniconda/Miniconda3-py311_23.5.2-0-Linux-x86_64.sh
#Follow-installation prompts
bash Miniconda3-py311_23.5.2-0-Linux-x86_64.sh
```

**Conda env creation**

```
conda install -n base -c conda-forge mamba
mamba env create -f envs/python_env.yaml
```

This will install dependencies into an isolated software environment, that has to be activated with

```
$ conda activate dual_annotation_env
```

## Using Python Virtual Environments

**Alternative option - python virtual environment**

```
pip install virtualenv
virtualenv dual_annotation_venv
python3 -m venv dual_annotation_venv
```

* Activate the environment and install dependencies

```
source dual_annotation_venv/bin/activate
pip install -r envs/requirements.txt
deactivate
```

This will install dependencies into an isolated software environment, that has to be activated with

```
$ source dual_annotation_venv/bin/activate
```

## Functionality

The core function, `compare_calls`, performs the following steps:

1. **Processing SMAP Files**: Read and parse `smap` files for dual and control annotations.
2. **Processing Aneuploidy Data**: Read aneuploidy data for both case and control.
3. **Processing CNV Data**: Stitch together CNV data, followed by reading it for both case and control.
4. **Intersecting Calls**: Processes and intersects CNVs and aneuploidies between case and control.
5. **Output**: Writes the results to various formats: Excel (`xlsx`), CSV (`csv`), and DOCX (`docx`).

## Parameters of `compare_calls`

- `dual_aneuploidy` (str): Path to the dual annotation aneuploidy file.
- `dual_smap` (str): Path to the dual annotation `smap` file.
- `dual_cnv` (str): Path to the dual annotation CNV file.
- `control_aneuploidy` (str): Path to the control aneuploidy file.
- `control_cnv` (str): Path to the control CNV file.
- `out_file` (str): Output file path where results will be saved.
- `case_id` (str): Identifier for the case sample.
- `control_id` (str): Identifier for the control sample.
- `celltype` (str): Type of cell being analyzed.
- `control_smap` (str): Path to the control `smap` file.
- `control_json` (str): JSON file path for the control data.
- `case_json` (str): JSON file path for the case data.
- `cnv_overlap_percentage` (float, default=0.3): Overlap percentage threshold for CNV processing.
- `aneuploidy_overlap_percentage` (float, default=0.5): Overlap percentage threshold for aneuploidy processing.
- `cnv_window` (int, default=1000): Window size for CNV processing.
- `cnv_stitch_window` (int, default=550000): Window size for stitching CNV data.


**Ensure that all filters are applied prior to downloading data**

### Cell QC - standard filtering recommendations​

| Filters                             | 400 Gbp                   | 1.5 Tbp                     | 5.0 Tbp                     |
|-------------------------------------|---------------------------|-----------------------------|-----------------------------|
| SV/CN confidence                     | Recommended               | Recommended                 | Recommended                 |
| SV Size (Insertions)                | 5kbp                      | 5kbp                        | 20kbp                       |
| SV Size (Deletions)                 | 7kbp                      | 7kbp                        | 50kbp                       |
| SV Size (Inversions)                | 50kbp                     | 70kbp                       | 100kbp                      |
| SV Size (Duplications)              | 50kbp                     | 140kbp                      | 100kbp                      |
| SV masking                          | Non-masked only           | Non-masked only             | Non-masked only             |
| CN Size                             | 10 Mb                     | 10 Mb                       | 10 Mb                       |
| CN Masking                          | All                       | All                         | All                         |
| Control db                          | 0% (absent from control db)| 0% (absent from control db) | 0% (absent from control db) |
| Found in paired control assembly    | No                        | No                          | No                          |
| Found in paired control molecules   | No                        | No                          | No                          |


Download the filtered CNV and Aneuploidy calls for the control/mock sample from Access. Download the smap, CNV, and Aneuploidy results from the Dual annotation results from Access. Place both results in a folder and unzip. To run the script activate an environment that contains the packages above and from the directory where `process_cnv_aneuploidy_calls.py` is located run the script.

Documentation describing data retrieval process from Access: [Genomic Data Comparison Tool](https://bionano.sharepoint.com/sites/clinical-scientific-affairs-global/_layouts/15/doc.aspx?sourcedoc={9bf3b32c-d1a9-4533-bd33-b09f399c71d7}&action=edit)

## RUNNING

Mock example where both the Mock and Dual annotation results were unzipped in the directory `dual_reporting`
```
python process_cnv_aneuploidy_calls.py --dual_aneuploidy data_package/GM16736_vs_Dilution/GM16736_-_Variant_Annotation_Pipeline_-_baseline_vs_0.1_dilution_7_13_2023_13_50_26_Aneuploidy.txt --dual_smap data_package/GM16736_vs_Dilution/GM16736_-_Variant_Annotation_Pipeline_-_baseline_vs_0.1_dilution_7_13_2023_13_50_26_Annotated_SV.smap --dual_cnv data_package/GM16736_vs_Dilution/GM16736_-_Variant_Annotation_Pipeline_-_baseline_vs_0.1_dilution_7_13_2023_13_50_26_CNV.txt --control_aneuploidy data_package/GM16736_vs_Dilution/RVP_dilution_VAF_0.1_seed_1_1.5tbp_7_13_2023_14_28_51_Aneuploidy.txt --control_cnv data_package/GM16736_vs_Dilution/RVP_dilution_VAF_0.1_seed_1_1.5tbp_7_13_2023_14_28_51_CNV.txt --control_smap data_package/GM16736_vs_Dilution/RVP_dilution_VAF_0.1_seed_1_1.5tbp_7_13_2023_14_28_51_Annotated_SV.smap --out_file GM16736_0.1_vs_baseline_test_update.xlsx --case_id 0.1_dilution --control_id GM16736_baseline --celltype HCM --control_json report_mock.json --case_json report_case.json 
```
