## Bionano - CNV and Aneuploidy dual annotation call integration
---
This script is designed to integrate and compare CNV and Aneuploidy calls from case and control samples with paired dual variant annotation results. This script will generate an .xlsx that contains
SVs, CNVs and Aneuploidies unique to the designated case.

## SETUP
---
Current implementation of this script utilizes the following packages
```
Python (3.10.4)
Pandas (1.4.3)
argparse (1.1)
numpy (1.23.2)
```

** Conda installation **

```
wget https://repo.anaconda.com/miniconda/Miniconda3-py39_4.12.0-Linux-x86_64.sh
#Follow-installation prompts
bash Miniconda3-py39_4.12.0-Linux-x86_64.sh
```

** Mamba installation **

```
conda install -n base -c conda-forge mamba
mamba env create -f envs/python_env.yaml
```

This will install dependencies into an isolated software environment, that has to be activated with

```
$ conda activate snakemake
$ snakemake --help
```

** Ensure that all filters are applied prior to downloading data **
Download the filtered smap, CNV and Aneuploidy calls for the control/mock sample from Access. Download the smap, CNV, and Aneuploidy results from the Dual annotation results from Access. For both the case and control sample assembly, select the informatics report and at the bottom of the report select "Download json". Place both results in a folder and unzip. To run the script activate an environment that contains the packages above and from the directory where `process_cnv_aneuploidy_calls.py` is located run the script.


## RUNNING

Mock example where both the Mock and Dual annotation results were unzipped in the directory `dual_reporting`
```
python process_cnv_aneuploidy_calls.py --dual_aneuploidy data_package/GM16736_vs_Dilution/GM16736_-_Variant_Annotation_Pipeline_-_baseline_vs_0.1_dilution_7_13_2023_13_50_26_Aneuploidy.txt --dual_smap data_package/GM16736_vs_Dilution/GM16736_-_Variant_Annotation_Pipeline_-_baseline_vs_0.1_dilution_7_13_2023_13_50_26_Annotated_SV.smap --dual_cnv data_package/GM16736_vs_Dilution/GM16736_-_Variant_Annotation_Pipeline_-_baseline_vs_0.1_dilution_7_13_2023_13_50_26_CNV.txt --control_aneuploidy data_package/GM16736_vs_Dilution/RVP_dilution_VAF_0.1_seed_1_1.5tbp_7_13_2023_14_28_51_Aneuploidy.txt --control_cnv data_package/GM16736_vs_Dilution/RVP_dilution_VAF_0.1_seed_1_1.5tbp_7_13_2023_14_28_51_CNV.txt --control_smap data_package/GM16736_vs_Dilution/RVP_dilution_VAF_0.1_seed_1_1.5tbp_7_13_2023_14_28_51_Annotated_SV.smap --out_file GM16736_0.1_vs_baseline_test_update.xlsx --case_id 0.1_dilution --control_id GM16736_baseline --celltype HCM --control_json report_mock.json --case_json report_case.json 
```
