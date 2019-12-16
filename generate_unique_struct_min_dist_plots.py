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
import VIU
from itertools import combinations, chain
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from seaborn import color_palette


LucksLabUtils_config.config("Quest_R2D2")

opts = OSU.getopts("", ["sample_sizes=", "outfile_pre=", "output_dir=", "reactivities_files=", "linker_seq=", "pol_fp=", "processors=", "MDS_processors="])
print opts
infiles = opts["--reactivities_files"].split(",")
outfile = opts["--outfile_pre"]
output_dir = OSU.create_directory(opts["--output_dir"])
reactivities_files = opts["--reactivities_files"].split(",")
linker_seq = opts["--linker_seq"]
sample_sizes = [int(s) for s in opts["--sample_sizes"].split(",")]
pol_fp = int(opts["--pol_fp"]) if "--pol_fp" in opts else 0
p = int(opts["--processors"]) if "--processors" in opts else 1
MDS_p = int(opts["--MDS_processors"]) if "--MDS_processors" in opts else 1

# setup counters, scaling functions, and output file header
unique_struct_nums = defaultdict(list)
unique_struct_dists = defaultdict(list)
noshape_struct_nums = defaultdict(list)
shape_struct_nums = defaultdict(list)
constrained_struct_nums = defaultdict(list)
scaling_fns = {"D": SU.invert_scale_rho_vec, "U": SU.scale_vec_avg1, "K": SU.cap_rho_or_ct_list}
with open(output_dir+"/"+outfile+".txt", "w") as f:
    f.write("Sample\t" + "\t".join([str(s) for s in sample_sizes]) + "\n")

# plotting setup
plt.style.use('seaborn-whitegrid')
fig = plt.figure(figsize = (3,3), dpi=300)

# iterate through reactivity files, general use case is only 1 reactivity file
for react_f in reactivities_files:
    file_name = re.findall("([^/]+).txt$", react_f)[0]
    output_file_prefix = output_dir + "/" + file_name
    pos, rho_full, theta, rho, seq, rc_flag, rc_sum, rc_untreated_sum, rc_treated_sum = SU.parse_reactivity_rho(react_f, linker_seq, output_file_prefix, endcut=-pol_fp)
    scaled_rhos = scaling_fns["K"](rho, 1)
    sampled_structs = set()

    # loop through sample sizes
    for s in range(len(sample_sizes)):
        e = sample_sizes[s] - sample_sizes[s-1] if s != 0 else sample_sizes[s]
        if e < 0:
            e = 0
        print "Currently sampling %s each"%(e)

        # Vanilla Sampling
        structs, structs_labels = SU.RNAstructure_sample(output_file_prefix, e, output_dir, label="noshape-%s"%(sample_sizes[s]*3), num_proc=1, wn_tag=file_name)
        sampled_structs.update(structs_labels)

        # Sampling with SHAPE constraints
        structs, structs_labels = SU.RNAstructure_sample(output_file_prefix, e, output_dir, shapefile=output_file_prefix+".rho", label="shape-%s"%(sample_sizes[s]*3), num_proc=1, wn_tag=file_name)
        sampled_structs.update(structs_labels)

        # Sampling with hard constraints
        XB = SU.get_indices_rho_gt_c(rho, 3.5, one_index=True)  # RNAstructure is 1-indexed
        SU.make_constraint_file(output_file_prefix+".con", [], XB, [], [], [], [], [], [])
        structs, structs_labels = SU.RNAstructure_sample(output_file_prefix, e, output_dir, constraintfile=output_file_prefix+".con", label="constrained_3.5-%s"%(sample_sizes[s]*3), num_proc=1, wn_tag=file_name)
        sampled_structs.update(structs_labels)
        sampled_structs_merged = SU.merge_labels(sampled_structs, to_string=False)

        # update output data structures
        unique_struct_nums[output_file_prefix].append(str(len(sampled_structs_merged)))
        # calculate counts
        label_counts = Counter([ss[1] for ss in sampled_structs_merged])
        label_counts_keys = label_counts.keys()
        reduced = [set([ik.split('-')[0] for ik in lk.split(',')]) for lk in label_counts_keys]
        noshape_only_keys = [rk[0] for rk in enumerate(reduced) if rk[1] == set(['noshape'])]
        noshape_only_count = sum([label_counts[noc[1]] for noc in enumerate(label_counts_keys) if noc[0] in noshape_only_keys])
        shape_only_keys = [rk[0] for rk in enumerate(reduced) if rk[1] == set(['shape'])]
        shape_only_count = sum([label_counts[soc[1]] for soc in enumerate(label_counts_keys) if soc[0] in shape_only_keys])
        constrained_only_keys = [rk[0] for rk in enumerate(reduced) if rk[1] == set(["constrained_"+str(3.5)])]
        constrained_only_count = sum([label_counts[coc[1]] for coc in enumerate(label_counts_keys) if coc[0] in constrained_only_keys])
        multiple_keys = [rk[0] for rk in enumerate(reduced) if len(rk[1]) > 1]
        multiple_count = sum([label_counts[mc[1]] for mc in enumerate(label_counts_keys) if mc[0] in multiple_keys])
        # saving counts of sampling method's unique structures
        noshape_struct_nums[output_file_prefix].append(str(noshape_only_count))
        shape_struct_nums[output_file_prefix].append(str(shape_only_count))
        constrained_struct_nums[output_file_prefix].append(str(constrained_only_count))

        # calculate distances
        sampled_structs_strings = [sl[0] for sl in sampled_structs_merged]
        ct_structs = [ss.split(",") for ss in sampled_structs_strings]
        SU.cts_to_file(ct_structs, seq, "%s/%s_%s_unique.ct"%(output_dir, outfile, sample_sizes[s]*3))
        binary_structs = SU.ct_struct_to_binary_vec(ct_structs)
        distances = []
        for bs in binary_structs:
            distances.append(SU.calc_bp_distance_vector_weighted(bs, scaled_rhos, scaling_func="K", invert_struct="D" != "K", paired_weight=0.8))
        unique_struct_dists[output_file_prefix].append(str(min(distances)))
    del sampled_structs, structs_labels, sampled_structs_strings, ct_structs

    # write outputfile
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
    del unique_struct_nums, unique_struct_dists, noshape_struct_nums, shape_struct_nums, constrained_struct_nums, 

    # clustering setup
    X = binary_structs
    Y = [",".join(sorted(set([al.split("-")[1] for al in a[1].split(",")]), key=int)) for a in sampled_structs_merged]
    Y_first_sampled = [yl.split(",")[0] for yl in Y]
    X_mats = SU.ct_struct_to_binary_mat([ssm[0].split(",") for ssm in sampled_structs_merged])
    first_color_dict = dict(zip(sorted(set(Y_first_sampled), key=int), color_palette("Spectral", len(sample_sizes))))
    del sampled_structs_merged

    with open("%s/%s_first_color_dict.txt" % (output_dir, outfile), "w") as f:
        f.write("\n".join("\t".join([k,str(v)]) for k,v in first_color_dict.items()) + "\n")
    with open("%s/%s_Y.txt" % (output_dir, outfile), "w") as f:
        f.write("\n".join(Y) + "\n")
    with open("%s/%s_Y_first_sampled.txt" % (output_dir, outfile), "w") as f:
        f.write("\n".join(Y_first_sampled) + "\n")

    # PCA
    principal_coords = SU.run_PCA(X, output_dir+"/"+outfile, center=False, scale_std=False)
    VIU.plot_PCA(principal_coords, Y_first_sampled, first_color_dict, output_dir+"/"+outfile, fig, plt, center=False, scale_std=False)

    # MDS
    print "Will use %s, %s processes" % (p, MDS_p)
    mds_mat_coords = SU.run_MDS_mat(X_mats, output_dir+"/"+outfile, p=p, MDS_p=MDS_p)
    VIU.plot_MDS(mds_mat_coords, Y_first_sampled, first_color_dict, output_dir+"/"+outfile, "mat", fig, plt)
    mds_ct_coords = SU.run_MDS_ct(X, output_dir+"/"+outfile, p=p, MDS_p=MDS_p)
    VIU.plot_MDS(mds_ct_coords, Y_first_sampled, first_color_dict, output_dir+"/"+outfile, "ct", fig, plt)

