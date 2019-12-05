"""
make_ct_to_PCA_MDS_plots.py
Takes in ct files, and creates PCA and MDS plots from them.

WARNING: Python's implementation of PCA may not give accurate coordinates for same structures
WARNING: MDS implementation in python may have same eigenvalue bug as cmdscale() in R
https://stat.ethz.ch/pipermail/r-sig-ecology/2010-July/001390.html

--input_cts		= Input ct file names with full path. Separated by ","
--input_names		= Input names to be used in each plot that describes the ct files
                          provided in "--intput_cts". Delimited by ","
--outfile		= Output file prefix
--output_dir		= Output directory

Author: Angela M Yu, 2014-2019
Version: 0.0.1

Copyright (C) 2019  Julius B. Lucks, Angela M Yu.
All rights reserved.
Distributed under the terms of the GNU General Public License, see 'LICENSE'.
"""
import OSU
import LucksLabUtils_config
import SU
import VIU
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from seaborn import color_palette
from itertools import repeat
from itertools import combinations, chain
import numpy

LucksLabUtils_config.config("Quest_R2D2")

opts = OSU.getopts("", ["input_cts=", "input_names=", "output_dir=", "outfile=", "processors="])
print opts
infiles = opts["--input_cts"].split(",")
input_names = opts["--input_names"].split(",")
output_dir = OSU.create_directory(opts["--output_dir"])
outfile = opts["--outfile"]
processors = int(opts["--processors"]) if "--processors" in opts else 1

output_file_prefix = output_dir + "/" + outfile
input_structs = set()
# load input ct files
for curr_ct_file, curr_name in zip(infiles, input_names):
    curr_structs = [(",".join(struct), name) for struct, name in zip(SU.get_ct_structs(curr_ct_file), repeat(curr_name))]
    input_structs.update(curr_structs)
unique_structs_merged = SU.merge_labels(input_structs, to_string=False)
del input_structs

# output unique structures and which dataset it was found
with open(output_file_prefix+".txt", "w") as f:
    f.write("\n".join("\t".join([k,v]) for k,v in unique_structs_merged) + "\n")

# set up combinations and colors for plotting
name_combinations = sorted([",".join([str(combo_e) for combo_e in combo]) for combo in chain(*map(lambda x: sorted(combinations(input_names, x)), range(1, len(input_names)+1)))])
color_dict = dict(zip(name_combinations, color_palette("Spectral", len(name_combinations))))
plt.style.use('seaborn-whitegrid')
fig = plt.figure(figsize = (3,3), dpi=300)

# set up X and Y
ct_structs = [ss[0].split(",") for ss in unique_structs_merged]
X = SU.ct_struct_to_binary_vec([ss[0].split(",") for ss in unique_structs_merged])
Y = [",".join(sorted(a[1].split(","))) for a in unique_structs_merged]  # sort labels in Y
X_mats = SU.ct_struct_to_binary_mat([ssm[0].split(",") for ssm in unique_structs_merged])
color_dict = dict(zip(sorted(set(Y)), color_palette("Spectral", len(name_combinations))))
del unique_structs_merged, name_combinations

with open("%s/%s_color_dict.txt" % (output_dir, outfile), "w") as f:
    f.write("\n".join("\t".join([k,str(v)]) for k,v in color_dict.items()) + "\n")
with open("%s/%s_Y.txt" % (output_dir, outfile), "w") as f:
    f.write("\n".join(Y) + "\n")

# PCA
principal_coords = SU.run_PCA(X, output_dir+"/"+outfile, center=False, scale_std=False)
VIU.plot_PCA(principal_coords, Y, color_dict, output_dir+"/"+outfile, fig, plt, center=False, scale_std=False)

# MDS
mds_ct_coords = SU.run_MDS_ct(X, output_dir+"/"+outfile, p=processors)
VIU.plot_MDS(mds_ct_coords, Y, color_dict, output_dir+"/"+outfile, "ct", fig, plt)
del X
mds_mat_coords = SU.run_MDS_mat(X_mats, output_dir+"/"+outfile, p=processors)
VIU.plot_MDS(mds_mat_coords, Y, color_dict, output_dir+"/"+outfile, "mat", fig, plt)

