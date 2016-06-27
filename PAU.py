"""
Panel Analysis Utilities (PAU)

Version: 0.0.0
Author: Angela M Yu
"""

import re
from collections import defaultdict
import SU
import NAU
import OSU
import cPickle as pickle
import operator
import subprocess
import time


def parse_input_panels(infiles, output_dir):
    """
    Parses input crystal ct files and experimental reactivities. This sets up the dictionary datastructure.
    It assumes different naming conventions to determine if the input is crystal structures or SHAPE-Seq data.
    """
    # infiles should be a directory with EITHER crystal structure ct files OR reactivities ct files
    parsed = {}
    for f in infiles:
        fname = re.findall("([^/]+_(\w+)_reactivities).txt$|([^/]+_(\w+)).ct$", f)  # Make sure to follow this regex
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


def train_constraint_model(crystal_dict, constrained_folds, cfr, react_rhos, rho_values, outdir, outname, weight=-1, scaling_func="D", cap_rhos=True):
    """
    Need to find the best rho value to sample with hard single-strand constraints
    Two-fold problem because:
    1. rhos are picked to find min distance by using as max_rho
    2. rhos to be picked to force as single stranded
    """
    if len(rho_values) == 1:
        output_dir = "/".join([outdir, "_".join([str(cfr), str(rho_values[0]), str(weight)])])
    else:
        output_dir = "/".join([outdir, "_".join([str(cfr), str(weight)])])

    scaling_fns = {"D": SU.invert_scale_rho_vec, "U": SU.scale_vec_avg1, "K": SU.cap_rho_or_ct_list}
    diffs_dict = dict.fromkeys(rho_values, 0)
    F_score_dict = dict.fromkeys(rho_values, 0)
    stats_dict = defaultdict(list)  # Python is weird sometimes
    min_dist_ind_dict = defaultdict(dict)
    for k in sorted(crystal_dict.keys()):
        for r in rho_values:
            # cap rho reactivities to a max value
            if cap_rhos:
                scaled_rhos = scaling_fns[scaling_func](react_rhos[k], r)
            else:
                scaled_rhos = scaling_fns[scaling_func](react_rhos[k])
            distances = []
            cf_dict_structs = [s.split(',') if isinstance(s, str) else s for s, l in constrained_folds[k]]

            # calculate distances between sampled structures and rho reactivities
            binary_structs = SU.ct_struct_to_binary_vec(cf_dict_structs)
            for s in binary_structs:
                distances.append(SU.calc_bp_distance_vector_weighted(s, scaled_rhos, scaling_func=scaling_func, invert_struct="U" in scaling_func, paired_weight=weight))
            min_distance = min(distances)
            min_dist_ind_dict[k][r] = [i for i, v in enumerate(distances) if v == min_distance]
            del binary_structs

            crystal_diffs = []
            F_score_list = []
            binary_struct_mats = SU.ct_struct_to_binary_mat(cf_dict_structs)
            for i in min_dist_ind_dict[k][r]:
                binary_struct_mat = binary_struct_mats[i]
                bm_stats = SU.calc_benchmark_statistics_matrix(binary_struct_mat, crystal_dict[k][2])
                crystal_diffs.append(SU.calc_bp_distance_matrix(binary_struct_mat, crystal_dict[k][2]))
                stats_dict[r].append(bm_stats)
                F_score_list.append(bm_stats["F_score"])
            avg_crystal_diffs = sum(crystal_diffs) / float(len(crystal_diffs))
            avg_F_score = sum(F_score_list) / float(len(F_score_list))
            diffs_dict[r] += avg_crystal_diffs
            F_score_dict[r] += avg_F_score
            with open(output_dir + "/" + outname + k + ".F", 'a') as f:
                f.write(str(k) + "\n")
                f.write("\t".join([str(r), str(avg_F_score)]) + "\n")
    output_train_model_stats(diffs_dict, stats_dict, min_dist_ind_dict, F_score_dict, crystal_dict, constrained_folds, output_dir)

    return diffs_dict, stats_dict, min_dist_ind_dict, F_score_dict


def load_train_model(rho_midpoint, constrain_val, paired_weight, reactivities, crystals, sample_n, react_rhos, structs_pickle_dir, output_dir, out_stat_dir, outname, scaling_func, cap_rhos, shape_slope, shape_intercept):
    """
    Runs the analysis of a (rho_midpoint, constrain_val, paired_weight) triple by loading in the sampled structures from sampling with hard constraints and constrain_val.
    """
    ret = []
    constrained_folds = {}
    for k in reactivities:
        if isinstance(reactivities[k][3], list):
            constrained_folds[k] = reactivities[k][3]
        else:
            constrained_folds[k] = []
    if isinstance(rho_midpoint, float):
        rho_values = [rho_midpoint]
    else:
        rho_values = rho_midpoint
    if len(rho_values) == 1:
        out_param_dir = "/".join([out_stat_dir, "_".join([str(constrain_val), str(rho_values[0]), str(paired_weight)])]) + "/"
    else:
        out_param_dir = "/".join([out_stat_dir, "_".join([str(constrain_val), str(paired_weight)])]) + "/"
    OSU.create_directory(out_param_dir)

    if constrain_val != "no_constrained":
        constrained_structs_dict = pickle.load(open("%sconstrained_folds_%s.p" % (structs_pickle_dir, str(constrain_val)), "rb"))
        for k in reactivities:
            constrained_structs = set(constrained_structs_dict[k])
            constrained_structs.update(set(constrained_folds[k]))
            constrained_structs = set(SU.merge_labels(list(constrained_structs), to_string=False))
            constrained_folds[k] = [(s.split(","), l) for s, l in constrained_structs]

    diffs_dict, stats_dict, min_dist_ind_dict, F_score_dict = train_constraint_model(crystals, constrained_folds, constrain_val, react_rhos, rho_values, out_stat_dir, outname + str(rho_midpoint), weight=paired_weight, scaling_func=scaling_func, cap_rhos=cap_rhos)

    best_F_score_dr = max(F_score_dict.iteritems(), key=operator.itemgetter(1))
    ret.append((rho_midpoint, constrain_val, paired_weight, best_F_score_dr[0], best_F_score_dr[1]))
    return ret


def load_train_model_helper(args):
    """
    Helper function for load_train_model
    """
    return load_train_model(*args)


def output_train_model_stats(all_diffs, all_diffs_stats, min_dist_struct_indices, F_scores, crystals, reactivities_structs, output_dir):
    """
    Output summary of the benchmarking
    """
    with open(output_dir + "/save_all_diffs.p", "wb") as f:
        pickle.dump(all_diffs, f)
    with open(output_dir + "/save_all_diffs_stats.p", "wb") as f:
        pickle.dump(all_diffs_stats, f)
    with open(output_dir + "/save_min_dist_struct_indices.p", "wb") as f:
        pickle.dump(min_dist_struct_indices, f)
    with open(output_dir + "/save_F_scores.p", "wb") as f:
        pickle.dump(F_scores, f)
    with open(output_dir + "/save_reactivities_structs.p", "wb") as f:
        pickle.dump(reactivities_structs, f)

    max_F_score = max(F_scores.values())
    max_F_keys = [k for k in sorted(F_scores) if F_scores[k] == max_F_score]

    for k in sorted(crystals.keys()):
        for mfk in max_F_keys:
            for mdi in min_dist_struct_indices[k][mfk]:
                outpre = output_dir+"/"+k+"_F"+str(mfk)+"_"+str(mdi)
                SU.ct_list_to_file(reactivities_structs[k][mdi][0], crystals[k][1], outpre+".ct")
                SU.runRNAstructure_CircleCompare(outpre+".ct", crystals[k][3], outpre+".ps")
                OSU.system_command("convert %s.ps %s.jpg" % (outpre, outpre))
    with open(output_dir + "/diffs_best_F.txt", 'w') as f:
        f.write("\t".join([str(m) for m in max_F_keys]) + "\n")
        f.write("Max avg F score: " + str(max_F_score) + "\n")
        for k in max_F_keys:
            f.write(str(k) + ": is the max_rho max_F_key\n")
            for ck in sorted(crystals.keys()):
                print "Crystal key: %s\n"%(str(ck))
                reactivities_labels = [reactivities_structs[ck][x][1] for x in min_dist_struct_indices[ck][k]]
                f.write("Methods: " + str(ck) + "\t" + str(reactivities_labels) + "\n")
            for stat_i in all_diffs_stats[k]:
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
