"""
sample_previous_lengths.py
Recalculates best structures when considering previous lengths' sampled structures

--seq_dir = Directory of sequence files
--out_dir = Output directory
--rhos_dir = Directory of rho reactivity files


Author: Angela M Yu, 2014-2018
Version: 0.0.1

Copyright (C) 2018 Julius B. Lucks, Angela M Yu.
All rights reserved.
Distributed under the terms of the GNU General Public License, see 'LICENSE'.
"""

import LucksLabUtils_config
import VIU
import OSU
import SU
import PCSU
import re
import glob
import pickle
import math

LucksLabUtils_config.config("Quest_R2D2")  # set up environment
OSU.system_command("echo $DATAPATH")
# parse command line arguments
opts = OSU.getopts("", ["100_repititions_prefix=", "out_dir=", "num_prev=", "scaling_func=", "weight_paired=", "rep_nums="])
print opts
repititions_prefix = opts["--100_repititions_prefix"]
out_dir = OSU.create_directory(opts["--out_dir"] + "/")
num_prev = int(opts["--num_prev"])
scaling_func = opts["--scaling_func"]
weight_paired = float(opts["--weight_paired"])
rep_nums = opts["--rep_nums"].split(",")

ct_dir = OSU.create_directory(out_dir + "ct_dir/")
draw_dir = OSU.create_directory(out_dir + "draw/")
efn2_dir = OSU.create_directory(out_dir + "efn2/")
pickle_dir = OSU.create_directory(out_dir + "pickles/")

# get reactivities and sampled structures from R2D2 iterations
rhos = {}
structs = {}
file_prefixes = {}
sequences = {}
for i in rep_nums: #range(1,101):
    curr_dir = "%s_%s"%(repititions_prefix, i)
    print curr_dir
    if i == rep_nums[0]:
        if not OSU.check_dir_exists("%s/rho_dir/"%(curr_dir)):
            OSU.system_command("tar -zxvf %s/results_except_draw.tgz -C %s ./rho_dir/"%(curr_dir, curr_dir))  # unzip rho reactivities only once
        for rf in glob.glob("%s/rho_dir/*.best_scaled_rho"%(curr_dir)):
            fname = re.findall("([^/]+).best_scaled_rho$", rf)[0]
            reactivities = SU.get_reactivity_from_file(rf)
            rhos[len(reactivities)] = reactivities
            file_prefixes[len(reactivities)] = fname
    if not OSU.check_dir_exists("%s/ct/"%(curr_dir)):
        OSU.system_command("tar -zxvf %s/results_except_draw.tgz -C %s ./ct/"%(curr_dir, curr_dir))  # unzip sampled structures at each length
    for ct in glob.glob("%s/ct/*_unique.ct"%(curr_dir)):
        ct_name = re.findall("([^/]+)_unique.ct$", ct)[0]
        curr_structs, curr_seqs = SU.get_ct_structs(ct, return_seq=True)
        structs[len(curr_structs[0])] = curr_structs
        sequences[len(curr_structs[0])] = "".join(curr_seqs[0])  # assumes all sequences the same per .ct file

# reads through .rho files found in rhos_dir
if set(rhos.keys()) != set(structs.keys()):
    raise Exception("reactivities and structure lengths are not the same")

results = {}
lengths = sorted(structs.keys())
print lengths
# TODO: parallelize
for li in range(len(lengths)):
    file_data_length_key = {}
    curr_structs = set()
    curr_length = lengths[li]
    for subset_lengths in lengths[max(li-num_prev, 0):(li+1)]:
        for s in structs[subset_lengths]:
            s_longer = s + ["0"] * (curr_length - len(s))
            curr_structs.add(tuple(s_longer))
    curr_structs = list(curr_structs)
    # calculate free energies
    SU.cts_to_file(curr_structs, sequences[curr_length], ct_dir+file_prefixes[curr_length]+"_unique.ct", shorten_name=True)
    SU.runRNAstructure_efn2(ct_dir+file_prefixes[curr_length]+"_unique.ct", efn2_dir + file_prefixes[curr_length] + ".efn2")
    free_energies = SU.get_free_energy_efn2(efn2_dir+ file_prefixes[curr_length] + ".efn2")
    file_data_length_key["free_energies"] = free_energies
    del free_energies
    # calculate distances
    binary_structs = SU.ct_struct_to_binary_vec(curr_structs)
    distances = []
    for s in binary_structs:
        distances.append(SU.calc_bp_distance_vector_weighted(s, rhos[curr_length], scaling_func=scaling_func, invert_struct="D" != scaling_func, paired_weight=weight_paired))
    file_data_length_key["distances"] = distances
    min_distance = min(distances)
    file_data_length_key["min_dist_indices"] = [i for i, v in enumerate(distances) if v == min_distance]
    file_data_length_key["min_distance"] = min_distance
    file_data_length_key["structs"] = curr_structs
    file_data_length_key["filename"] = file_prefixes[curr_length]
    file_data_length_key["rho"] = rhos[curr_length]
    file_data_length_key["rc_flag"] = 1  # TODO: make this work
    # save data
    pickle.dump(file_data_length_key, open("%sfile_data_%s.p"%(pickle_dir, curr_length), "wb"))
    results[curr_length] = file_data_length_key

PCSU.generate_DG_output(results, output_dir=out_dir)

zero_padding = int(math.floor(math.log10(lengths[-1])) + 1)
varna_num = range(1, len(results)+1)
for length, vn in zip(lengths, varna_num):
    PCSU.generate_best_struct_images(results[length], length, lengths[-1], vn, zero_padding, draw_dir, ct_dir, draw_all=True, most_count_tie_break=False)


