"""
Panel Analysis Utilities (PAU)

Version: 0.0.1
Author: Angela M Yu

Copyright (C) 2016  Julius B. Lucks, Angela M Yu.
All rights reserved.
Distributed under the terms of the GNU General Public License, see 'LICENSE'.
"""

import re
from collections import defaultdict
import SU
import NAU
import OSU
import cPickle as pickle
import subprocess
import time


def parse_input_panels(infiles, output_dir):
    """
    Parses input crystal ct files and experimental reactivities (infiles). This sets up the dictionary datastructure.
    It assumes distinct naming conventions to determine if the input is crystal structures or SHAPE-Seq data.
    Also requires output_dir
    """
    # infiles should be a directory with EITHER crystal structure ct files OR reactivities ct files
    parsed = {}
    for f in infiles:
        fname = re.findall('([^/]+_([\\w]+[-]?[\\w]*)_reactivities).txt$|([^/]+_(\\w+)).ct$', f)  # Make sure to follow this regex
        if fname[0][1] != '':  # reactivities file
            output_file_prefix = output_dir + "/" + fname[0][0]
            pos, rho_orig, theta, rho_cut, seq, rc_flag, rc_sum, ut_sum, t_sum = SU.parse_reactivity_rho(f, "", output_file_prefix)
            parsed[fname[0][1]] = [rho_orig, seq, output_file_prefix]
        elif fname[0][3] != '':  # crystal structure
            struct = SU.get_ct_structs(f)
            # rho here is scaled and in matrix format (not upper triangular)
            rho = SU.ct_struct_to_binary_vec(struct)
            seq = NAU.get_seq_from_ct(f)
            struct_mat = SU.ct_struct_to_binary_mat(struct)
            parsed[fname[0][3]] = [rho[0], seq[0], struct_mat[0], f, struct[0]]
    return parsed


def train_constraint_model(crystal_dict, constrained_folds, cfr, react_rhos, rho_value, outdir, outname, weight=-1, scaling_func="D", cap_rhos=True):
    """
    Need to find the best rho value to sample with hard single-strand constraints
    Two-fold problem because:
    1. rhos are picked to find min distance by using as max_rho
    2. rhos to be picked to force as single stranded
    """
    output_dir = "/".join([outdir, "_".join([str(cfr), str(rho_value), str(weight)])])
    scaling_fns = {"D": SU.invert_scale_rho_vec, "U": SU.scale_vec_avg1, "K": SU.cap_rho_or_ct_list}
    F_score = 0
    all_stats = defaultdict(list)
    min_dist_ind_dict = defaultdict(list)
    # Obsolete now that we are calling this function for every parameter set
    for k in sorted(react_rhos.keys()):
        ck = k.split('-')[0]  # corresponding crystal key from the reactivity key
        # cap rho reactivities to a max value
        if cap_rhos:
            scaled_rhos = scaling_fns[scaling_func](react_rhos[k], rho_value)
        else:
            scaled_rhos = scaling_fns[scaling_func](react_rhos[k])
        distances = []
        cf_structs = [s.split(',') if isinstance(s, str) else s for s, l in constrained_folds[k]]

        # calculate distances between sampled structures and rho reactivities
        binary_structs = SU.ct_struct_to_binary_vec(cf_structs)
        for s in binary_structs:
            distances.append(SU.calc_bp_distance_vector_weighted(s, scaled_rhos, scaling_func=scaling_func, invert_struct="D" != scaling_func, paired_weight=weight))
        min_distance = min(distances)
        min_dist_ind_dict[k] = [i for i, v in enumerate(distances) if v == min_distance]
        del binary_structs

        # calculate and output results for this parameter set
        crystal_diffs = []
        F_score_list = []
        for i in min_dist_ind_dict[k]:
            print "\nStruct key %s min distance index %s" % (k, i)
            binary_struct_mat = SU.ct_struct_to_binary_mat(cf_structs[i])
            bm_stats = SU.calc_benchmark_statistics_matrix(binary_struct_mat, crystal_dict[ck][2])
            crystal_diffs.append(SU.calc_bp_distance_matrix(binary_struct_mat, crystal_dict[ck][2]))
            all_stats[k].append(bm_stats)
            F_score_list.append(bm_stats["F_score"])
        avg_F_score = sum(F_score_list) / float(len(F_score_list))
        F_score += avg_F_score
    output_train_model_stats(all_stats, min_dist_ind_dict, F_score, crystal_dict, constrained_folds, output_dir)

    return all_stats, min_dist_ind_dict, F_score


def load_train_model(rho_midpoint, constrain_val, paired_weight, reactivities, crystals, sample_n, react_rhos, structs_pickle_dir, output_dir, out_stat_dir, outname, scaling_func, cap_rhos, shape_slope, shape_intercept):
    """
    Runs the analysis of a (rho_midpoint, constrain_val, paired_weight) triple by loading in the sampled structures from sampling with hard constraints and constrain_val.
    """
    constrained_folds = {}
    # Obselete since we use 1 reactivity at a time
    for k in reactivities:
        if isinstance(reactivities[k][3], list):
            constrained_folds[k] = reactivities[k][3]
        else:
            constrained_folds[k] = []
    out_param_dir = "/".join([out_stat_dir, "_".join([str(constrain_val), str(rho_midpoint), str(paired_weight)])]) + "/"
    OSU.create_directory(out_param_dir)

    if constrain_val != "no_constrained":
        # load sampled folds with hard constrain c
        constrained_structs_dict = pickle.load(open("%sconstrained_folds_%s.p" % (structs_pickle_dir, str(constrain_val)), "rb"))
        # Obselete since we use 1 reactivity at a time
        for k in reactivities:
            constrained_structs = set(constrained_structs_dict[k])
            constrained_structs.update(set(constrained_folds[k]))
            constrained_structs = set(SU.merge_labels(list(constrained_structs), to_string=False))
            constrained_folds[k] = [(s.split(","), l) for s, l in constrained_structs]

    # call training routine
    stats_dict, min_dist_ind_dict, F_score = train_constraint_model(crystals, constrained_folds, constrain_val, react_rhos, rho_midpoint, out_stat_dir, outname + str(rho_midpoint), weight=paired_weight, scaling_func=scaling_func, cap_rhos=cap_rhos)

    return [[rho_midpoint, constrain_val, paired_weight, F_score]]


def load_train_model_helper(args):
    """
    Helper function for load_train_model
    """
    return load_train_model(*args)


def output_train_model_stats(all_diffs_stats, min_dist_struct_indices, F_score, crystals, reactivities_structs, output_dir):
    """
    Output summary of the benchmarking
    """
    # Pickle benchmark statistics and indices of structures with the minimum distance
    with open(output_dir + "/save_all_diffs_stats.p", "wb") as f:
        pickle.dump(all_diffs_stats, f)
    with open(output_dir + "/save_min_dist_struct_indices.p", "wb") as f:
        pickle.dump(min_dist_struct_indices, f)

    # Generate circle compare diagrams of minimum distance structures
    for k in sorted(reactivities_structs.keys()):
        ck = k.split('-')[0]  # corresponding crystal key from the reactivity key
        for mdi in min_dist_struct_indices[k]:
            outpre = output_dir+"/"+k+str(mdi)
            SU.ct_list_to_file(reactivities_structs[k][mdi][0], crystals[ck][1], outpre+".ct")
            SU.runRNAstructure_CircleCompare(outpre+".ct", crystals[ck][3], outpre+".ps")
            OSU.system_command("convert %s.ps %s.jpg" % (outpre, outpre))

    # Write out information on highest F score structures
    with open(output_dir + "/diffs_best_F.txt", 'w') as f:
        f.write("Max avg F score: " + str(F_score) + "\n")
        for k in sorted(reactivities_structs.keys()):
            reactivities_labels = [reactivities_structs[k][x][1] for x in min_dist_struct_indices[k]]
            f.write("Methods: " + str(k) + "\t" + str(reactivities_labels) + "\n")
        for k, stat_i in all_diffs_stats.items():
            f.write("Structure key: %s\n" % (k))
            f.write(str(stat_i) + "\n")
    return


def wait_for_jcoll(jcoll, out_dir, wait=30):
    """
    Makes current process wait until no jobs are left in a job collection jcoll
    """
    wait_jlist = True
    while wait_jlist:
        jlist_jcoll = subprocess.Popen(['/opt/voyager/nbs/bin/jlist', '-jcoll', jcoll], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        out, err = jlist_jcoll.communicate()
        if "listempty" in out:
            wait_jlist = False
            break
        OSU.system_command("echo \"waiting\n%s\n%s\n\" >> %sjcoll_waiting.txt" % (out, err, out_dir))
        time.sleep(wait)
    OSU.system_command("echo \"%s\n\" >> %sjcoll_waiting.txt" % (time.localtime(), out_dir))


def wait_jcoll_finish_any(jcoll, out_dir, max_jobs, wait=30):
    """
    Waits until any job in the job collection is done before returning the number of available jobs
    """
    wait_jlist = True
    while wait_jlist:
        num_running = count_jcoll_remaining(jcoll, out_dir)
        if num_running < max_jobs:
            wait_jlist = False
            OSU.system_command("echo \"jobs available\n%s\n%s\n\" >> %sjcoll_waiting.txt" % (num_running, max_jobs, out_dir))
            break
        OSU.system_command("echo \"waiting jobs available\n%s\n%s\n\" >> %sjcoll_waiting.txt" % (num_running, max_jobs, out_dir))
        time.sleep(wait)
    OSU.system_command("echo \"%s\n\" >> %sjcoll_waiting.txt" % (time.localtime(), out_dir))
    return max_jobs - num_running


def count_jcoll_remaining(jcoll, out_dir):
    """
    Counts number of jobs remaining in job collection jcoll
    """
    jlist_jcoll = subprocess.Popen(['/opt/voyager/nbs/bin/jlist', '-jcoll', jcoll], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = jlist_jcoll.communicate()
    count_jobs = 0 if "listempty" in out else out.count("\n") - 5
    return count_jobs


def check_job_on_queue(job_name):
    """
    Returns boolean of if job_name is on the queue
    """
    jlist = subprocess.Popen(['/opt/voyager/nbs/bin/jlist'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = jlist.communicate()
    return job_name in out
