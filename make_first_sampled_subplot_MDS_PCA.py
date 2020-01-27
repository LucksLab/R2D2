"""
make_first_sampled_subplot_MDS_PCA.py
Generates plots of previously computed MDS and PCA coordinates, but only plotting subsections of the data for better visualization.
To be used after generate_unique_struct_min_dist_plots.py

--Y_first_sampled	= First sampled label file
--MDS_ct_coords		= MDS ct coordinates file
--MDS_mat_coords        = MDS mat coordinates file
--PCA_coords	        = PCA ct coordinates file
--outfile_pre		= Output file prefix
--output_dir		= Output directory

Author: Angela M Yu, 2014-2019
Version: 0.0.1

Copyright (C) 2019  Julius B. Lucks, Angela M Yu.
All rights reserved.
Distributed under the terms of the GNU General Public License, see 'LICENSE'.
"""
import OSU
import LucksLabUtils_config
import VIU
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from seaborn import color_palette
import numpy as np


LucksLabUtils_config.config("Quest_R2D2")

opts = OSU.getopts("", ["Y_first_sampled=", "MDS_ct_coords=", "MDS_mat_coords=", "PCA_coords=", "outfile_pre=", "output_dir=", "first_color_dict="])
print opts
filename_Y_first_sampled = opts["--Y_first_sampled"]
filename_MDS_ct_coords = opts["--MDS_ct_coords"]
filename_MDS_mat_coords = opts["--MDS_mat_coords"]
filename_PCA_coords = opts["--PCA_coords"]
outfile = opts["--outfile_pre"]
output_dir = OSU.create_directory(opts["--output_dir"])
filename_first_color_dict = opts["--first_color_dict"]

# plotting setup
plt.style.use('seaborn-whitegrid')
fig = plt.figure(figsize = (3,3), dpi=300)

# load input
with open(filename_Y_first_sampled, "r") as f:
    Y_first_sampled = [int(val.rstrip()) for val in f.readlines()]

with open(filename_MDS_ct_coords, "r") as f:
    MDS_ct_coords = [[float(ele) for ele in line.split()] for line in f.readlines()]

with open(filename_MDS_mat_coords, "r") as f:
    MDS_mat_coords = [[float(ele) for ele in line.split()] for line in f.readlines()]

with open(filename_PCA_coords, "r") as f:
    PCA_coords = [[float(ele) for ele in line.split()] for line in f.readlines()]

with open(filename_first_color_dict, "r") as f:
    first_color_dict = {int(k): tuple(float(cv) for cv in col[1:-2].split(", ")) for k, col in [line.split("\t") for line in f.readlines()]}

# min and max coordinate values
xlim_PCA = (min([x for x, y in PCA_coords]) - 1, max([x for x, y in PCA_coords]) + 1)
ylim_PCA = (min([y for x, y in PCA_coords]) - 1, max([y for x, y in PCA_coords]) + 1)

xlim_MDS_ct = (min([x for x, y in MDS_ct_coords]) - 1, max([x for x, y in MDS_ct_coords]) + 1)
ylim_MDS_ct = (min([y for x, y in MDS_ct_coords]) - 1, max([y for x, y in MDS_ct_coords]) + 1)

xlim_MDS_mat = (min([x for x, y in MDS_mat_coords]) - 1, max([x for x, y in MDS_mat_coords]) + 1)
ylim_MDS_mat = (min([y for x, y in MDS_mat_coords]) - 1, max([y for x, y in MDS_mat_coords]) + 1)

# iterate through increasing keys
ordered_keys = sorted(first_color_dict.keys())
for i in xrange(1, len(ordered_keys)+1):
    # subselect data
    sub_keys = ordered_keys[:i]
    MDS_mat_coords_sub, MDS_ct_coords_sub, Y_first_sampled_sub, PCA_coords_sub = zip(*[ele for ele in zip(MDS_mat_coords, MDS_ct_coords, Y_first_sampled, PCA_coords) if ele[2] in sub_keys])

    # plot
    VIU.plot_MDS(np.matrix(MDS_mat_coords_sub), Y_first_sampled_sub, first_color_dict, "%s/%s_%s" % (output_dir, outfile, "_".join([str(sk) for sk in sub_keys])), "mat", fig, plt, xlim=xlim_MDS_mat, ylim=ylim_MDS_mat)
    VIU.plot_MDS(np.matrix(MDS_ct_coords_sub), Y_first_sampled_sub, first_color_dict, "%s/%s_%s" % (output_dir, outfile, "_".join([str(sk) for sk in sub_keys])), "ct", fig, plt, xlim=xlim_MDS_ct, ylim=ylim_MDS_ct)
    VIU.plot_PCA(np.matrix(PCA_coords_sub), Y_first_sampled_sub, first_color_dict, "%s/%s_%s" % (output_dir, outfile, "_".join([str(sk) for sk in sub_keys])), fig, plt, center=False, scale_std=False, xlim=xlim_PCA, ylim=ylim_PCA)

