import os
import numpy as np
from pathlib import Path
import pandas as pd


work_dir = config["work_dir"]
log_dir = "Log"
symlink_dir = "Symlink"
alignment_dir = "Alignment"
maf_dir = "MAF"
insert_dir = "Insert_Size"
tmp_blat_dir = "Tmp_Blat"
blat_final_dir = "Blat"
ref_data_dir = "Reference_Data"
filtered_dir = "Filtered"

bam_source_dir = config["bam_source_dir"]
maf_source_dir = config["maf_source_dir"]
insert_size_source_dir = config["insert_size_source_dir"]
sample_list = config["sample_list"]
scatter_count = config["scatter_count"]


genome=config["genome"]
blat_binary = config["blat_binary"]
faToTwoBit_binary = config["faToTwoBit_binary"]
blat_filter_script_dir = config["blat_filter_script_dir"]
blat_filter_container = config["blat_filter_container"]
