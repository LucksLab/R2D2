"""
generate_unique_structures.py
Generates plots of number of unique structures and minimum distance of the sampled structures given reactivities files.

--reactivities_files	= Reactivities files to use to direct sampling and minimum distance structure(s)
--outfile_pre		= Output file prefix
--output_dir		= Output directory
--linker_seq		= Linker sequence
--sample_sizes		= Comma delimited list of sample sizes

Author: Angela M Yu, 2014-2017
Version: 0.0.1

Copyright (C) 2016, 2017  Julius B. Lucks, Angela M Yu.
All rights reserved.
Distributed under the terms of the GNU General Public License, see 'LICENSE'.
"""
import OSU
import LucksLabUtils_config
from collections import defaultdict, Counter
import re
import SU
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


LucksLabUtils_config.config("Quest_R2D2")

opts = OSU.getopts("", ["sample_sizes=", "outfile_pre=", "output_dir=", "reactivities_files=", "linker_seq="])

infiles = opts["--reactivities_files"].split(",")
outfile = opts["--outfile_pre"]
output_dir = OSU.create_directory(opts["--output_dir"])
reactivities_files = opts["--reactivities_files"].split(",")
linker_seq = opts["--linker_seq"]
sample_sizes = [int(s) for s in opts["--sample_sizes"].split(",")]

unique_struct_nums = defaultdict(list)
unique_struct_dists = defaultdict(list)
noshape_struct_nums = defaultdict(list)
shape_struct_nums = defaultdict(list)
constrained_struct_nums = defaultdict(list)
scaling_fns = {"D": SU.invert_scale_rho_vec, "U": SU.scale_vec_avg1, "K": SU.cap_rho_or_ct_list}
with open(output_dir+"/"+outfile+".txt", "w") as f:
    f.write("Sample\t" + "\t".join([str(s) for s in sample_sizes]) + "\n")

for react_f in reactivities_files:
    file_name = re.findall("([^/]+).txt$", react_f)[0]
    output_file_prefix = output_dir + "/" + file_name
    pos, rho_full, theta, rho, seq, rc_flag, rc_sum, rc_untreated_sum, rc_treated_sum = SU.parse_reactivity_rho(react_f, linker_seq, output_file_prefix)
    scaled_rhos = scaling_fns["K"](rho, 1)
    sampled_structs = set()

    for s in range(len(sample_sizes)):
        e = sample_sizes[s] - sample_sizes[s-1] if s != 0 else sample_sizes[s]
        if e < 0:
            e = 0
        print "Currently sampling %s each"%(e)

        # Vanilla Sampling
        structs, structs_labels = SU.RNAstructure_sample(output_file_prefix, e, output_dir, label="noshape", num_proc=1)
        sampled_structs.update(structs_labels)

        # Sampling with SHAPE constraints
        structs, structs_labels = SU.RNAstructure_sample(output_file_prefix, e, output_dir, shapefile=output_file_prefix+".rho", label="shape", num_proc=1)
        sampled_structs.update(structs_labels)

        # Sampling with hard constraints
        XB = SU.get_indices_rho_gt_c(rho, 3.5, one_index=True)  # RNAstructure is 1-indexed
        SU.make_constraint_file(output_file_prefix+".con", [], XB, [], [], [], [], [], [])
        structs, structs_labels = SU.RNAstructure_sample(output_file_prefix, e, output_dir, constraintfile=output_file_prefix+".con", label="constrained_"+str(3.5), num_proc=1)
        sampled_structs.update(structs_labels)
        sampled_structs_merged = SU.merge_labels(sampled_structs, to_string=False)

        # update output data structures
        unique_struct_nums[output_file_prefix].append(str(len(sampled_structs_merged)))
        label_counts = Counter([ss[1] for ss in sampled_structs_merged])
        noshape_struct_nums[output_file_prefix].append(str(label_counts["noshape"]))
        shape_struct_nums[output_file_prefix].append(str(label_counts["shape"]))
        constrained_struct_nums[output_file_prefix].append(str(label_counts["constrained_"+str(3.5)]))

        sampled_structs_strings = [sl[0] for sl in sampled_structs_merged]
        binary_structs = SU.ct_struct_to_binary_vec([ss.split(",") for ss in sampled_structs_strings])
        distances = []
        for s in binary_structs:
            distances.append(SU.calc_bp_distance_vector_weighted(s, scaled_rhos, scaling_func="K", invert_struct="D" != "K", paired_weight=0.8))

        unique_struct_dists[output_file_prefix].append(str(min(distances)))

    with open(output_dir+"/"+outfile+".txt", "a") as f:
        f.write(file_name + "_total_unique_structs\t")
        f.write("\t".join(unique_struct_nums[output_file_prefix]) + "\n")
        f.write(file_name + "_min_dist\t")
        f.write("\t".join(unique_struct_dists[output_file_prefix]) + "\n")
        f.write(file_name + "_noshape_unique_structs\t")
        f.write("\t".join(noshape_struct_nums[output_file_prefix]) + "\n")
        f.write(file_name + "_shape_unique_structs\t")
        f.write("\t".join(shape_struct_nums[output_file_prefix]) + "\n")
        f.write(file_name + "_constrained_unique_structs\t")
        f.write("\t".join(constrained_struct_nums[output_file_prefix]) + "\n")

    plt.plot([x*3 for x in sample_sizes], unique_struct_nums[output_file_prefix]) #, marker='o')
    plt.xlabel("Sample size")
    plt.ylabel("Number of unique structures")
    plt.title("Sampling RNA structures")
    plt.savefig(output_dir+"/"+outfile+'.png')

    plt.gcf().clear()

    plt.plot([x*3 for x in sample_sizes], unique_struct_dists[output_file_prefix]) #, marker='o')
    plt.xlabel("Sample size")
    plt.ylabel("Min distance")
    plt.title("Sampling RNA structures")
    plt.savefig(output_dir+"/"+outfile+'_distances.png')
