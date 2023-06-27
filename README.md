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
```

** Ensure that all filters are applied prior to downloading data **
Download the filtered CNV and Aneuploidy calls for the control/mock sample from Access. Download the smap, CNV, and Aneuploidy results from the Dual annotation results from Access. Place both results in a folder and unzip. To run the script activate an environment that contains the packages above and from the directory where `process_cnv_aneuploidy_calls.py` is located run the script.


## RUNNING

Mock example where both the Mock and Dual annotation results were unzipped in the directory `dual_reporting`
```
python process_cnv_aneuploidy_calls.py --dual_aneuploidy dual_reporting/HsMM-059-D2/HsMM-059-D2-RNP-versus-Mock_6_22_2023_21_23_22/HsMM-059-D2-RNP-versus-Mock_6_22_2023_21_23_22_Aneuploidy.txt --dual_smap dual_reporting/HsMM-059-D2/HsMM-059-D2-RNP-versus-Mock_6_22_2023_21_23_22/HsMM-059-D2-RNP-versus-Mock_6_22_2023_21_23_22_Annotated_SV.smap --dual_cnv dual_reporting/HsMM-059-D2/HsMM-059-D2-RNP-versus-Mock_6_22_2023_21_23_22/HsMM-059-D2-RNP-versus-Mock_6_22_2023_21_23_22_CNV.txt --control_aneuploidy dual_reporting/HsMM-059-D2/HsMM-059-Mock_RVP_6_22_2023_21_26_4/HsMM-059-Mock_RVP_6_22_2023_21_26_4_Aneuploidy.txt --control_cnv dual_reporting/HsMM-059-D2/HsMM-059-Mock_RVP_6_22_2023_21_26_4/HsMM-059-Mock_RVP_6_22_2023_21_26_4_CNV.txt --out_file test_HsMM_059.xlsx
```