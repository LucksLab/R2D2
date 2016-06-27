"""
Benchmarking code R2D2 (Reconstructing RNA Dynamics from Data) using grid
search on different parameter sets.
Version: 0.0.0
Author: Angela M Yu, 2014-2016

There may be different naming conventions of the input files, so the code
assumes input files to be of the form:
 *_(\w+).ct and *_(\w+)_reactivities.txt
Matching (\w+) will be compared and (\w+) is assumed to be a unique tag.
"""

import glob
import SU
import OSU
import NAU
import re
import LucksLabUtils_config
import PAU
from itertools import cycle
from collections import namedtuple
import cPickle as pickle

# setup environment variables specific to the ICSE cluster at Cornell
LucksLabUtils_config.config("ICSE")
opts = OSU.getopts("o:c:r:n:p:",
                   ["noshape", "shape", "constrain", "pseudoknots", "keys=",
                    "best_dist_func=", "scaling_func=", "cluster_flag=",
                    "job_name=", "sub_proc=", "arg_slice=", "load_results=",
                    "generate_structs=", "structs_pickle_dir=", "cap_rhos=",
                    "add_crystal_struct", "shape_intercept=", "shape_slope=",
                    "restart"])
print opts

reactivities_files = glob.glob(opts['-r'])
crystal_files = glob.glob(opts['-c'])
output_dir = opts['-o']
sample_n = int(opts['-n'])
num_proc = int(opts['-p'])
best_dist_func = opts["--best_dist_func"]
scaling_func = opts["--scaling_func"]
cluster_flag = opts["--cluster_flag"] == "True"
job_name = opts["--job_name"]
sub_proc = opts["--sub_proc"] == "True"
load_results = opts["--load_results"] == "True"
generate_structs = opts["--generate_structs"] == "True"
structs_pickle_dir = opts["--structs_pickle_dir"]
cap_rhos = opts["--cap_rhos"] == "True"
shape_intercept = float(opts["--shape_intercept"]) \
                    if "--shape_intercept" in opts else -0.3
shape_slope = float(opts["--shape_slope"]) if "--shape_slope" in opts else 1.1

restrict = []
if '--keys' in opts:  # remove these keys from analysis
    restrict = opts['--keys'].split(",")

OSU.create_directory(output_dir)
reactivities = PAU.parse_input_panels(reactivities_files, output_dir)
crystals = PAU.parse_input_panels(crystal_files, output_dir)

for k in restrict:
    crystals.pop(k, None)
    reactivities.pop(k, None)

if set(crystals.keys()) != set(reactivities.keys()):
    raise Exception("Keys of crystal structures and reactivities not equal")

outname = "diffs_all.txt"

# Parallelize
# also strips seq and updates rhos accordingly
for k, v in reactivities.iteritems():
    react_rhos = v[0]
    react_seq = v[1]
    cryst_seq = crystals[k][1]
    ind_of_match = react_seq.index(cryst_seq)
    end_of_match = ind_of_match + len(cryst_seq)
    # renormalize rhos based on matching regions
    react_rhos = SU.recalc_rhos(react_rhos, ind_of_match, end_of_match)

    # remove 5' from rhos used in further in the analysis
    # reactivities[k][0] = react_rhos[ind_of_match:]
    reactivities[k][0] = react_rhos
    reactivities[k][1] = react_seq
    SU.rhos_list_to_file(reactivities[k][0], reactivities[k][2]+".rho")

    sampling_opts_string = []
    if generate_structs:
        seqfile = NAU.make_seq(cryst_seq, output_dir+"temp.seq")
        sampled_structs = set()
        if "--noshape" in opts:
            structs, structs_labels = SU.RNAstructure_sample(reactivities[k][2], sample_n, output_dir, seqfile, shapefile="", constraintfile="", label="noshape", num_proc=num_proc, shape_slope=shape_slope, shape_intercept=shape_intercept)
            sampled_structs.update(structs_labels)
            print "after noshape: " + str(len(sampled_structs))
            sampling_opts_string.append("--noshape")
        if "--shape" in opts:
            structs, structs_labels = SU.RNAstructure_sample(reactivities[k][2], sample_n, output_dir, seqfile, shapefile=reactivities[k][2]+".rho", constraintfile="", label="shape", num_proc=num_proc, shape_slope=shape_slope, shape_intercept=shape_intercept)
            sampled_structs.update(structs_labels)
            print "after shape: " + str(len(sampled_structs))
            sampling_opts_string.append("--shape")
        sampled_structs = SU.merge_labels(list(sampled_structs), to_string=False)
        print "number of sampled structs: " + str(len(sampled_structs))
        reactivities[k].append(sampled_structs)
        OSU.remove_file(output_dir+"temp.seq")

out_stat_dir = output_dir + "/stats_out/"
OSU.create_directory(out_stat_dir)
OSU.create_directory(structs_pickle_dir)
CRW_pair = namedtuple("CRW_pair", ["c", "mr", "w"])

react_rhos = dict((k, reactivities[k][0]) for k in reactivities)
constrain_rho_F = {}

# Reduce search space
rho_midpoints = [0.1 * i for i in range(7, 61)]
weights = [0.1 * i for i in range(11)]
print "Parameter values to test: "
print "rho_max and rho_c: " + str(rho_midpoints)
print "weights: " + str(weights)

if '--constrain' not in opts:
    constrain_vals = ["no_constrained"]
    out_sub_dir = out_stat_dir + "no_constrained/"
elif len(opts['--constrain']) > 0:
    constrained_vals = [float(c) for c in opts['--constrain'].split(',')]
    out_sub_dir = out_stat_dir + "constrained/"
    sampling_opts_string.append("--constrain %s" % (opts['--constrain']))
else:
    constrain_vals = rho_midpoints
    out_sub_dir = out_stat_dir + "rho_midpoint_constrained/"
    OSU.create_directory(out_sub_dir)
    sampling_opts_string.append("--constrain")
sampling_opts_string = " ".join(sampling_opts_string)

# Only generate a set of structures once to be used in the benchmarking
if generate_structs:
    pickle.dump(crystals, open(structs_pickle_dir + "/save_crystals.p", "wb"))
    pickle.dump(reactivities, open(structs_pickle_dir + "/save_reactivities.p", "wb"))
    for constrain_val in constrain_vals:
        constrained_folds = {}
        for k in reactivities:
            XB = SU.get_indices_rho_gt_c(reactivities[k][0], constrain_val, one_index=True)
            SU.make_constraint_file(structs_pickle_dir + "/temp_c%s.con" % (str(constrain_val)), [], XB, [], [], [], [], [], [])
            seqfile = NAU.make_seq(crystals[k][1], structs_pickle_dir+"temp_c%s.seq" % (str(constrain_val)))
            sampled_counts, sampled = SU.RNAstructure_sample(reactivities[k][2], sample_n, structs_pickle_dir, seqfile, shapefile="", constraintfile=structs_pickle_dir + "/temp_c%s.con" % (str(constrain_val)), label="constrain_"+str(constrain_val), num_proc=1, shape_slope=shape_slope, shape_intercept=shape_intercept)
            constrained_folds[k] = sampled
            OSU.remove_files([structs_pickle_dir + "/temp_c%s.con" % (str(constrain_val)), structs_pickle_dir+"temp_c%s.seq" % (str(constrain_val))])
        pickle.dump(constrained_folds, open("%sconstrained_folds_%s.p" % (structs_pickle_dir, str(constrain_val)), "wb"))
else:
    with open(structs_pickle_dir + "/save_crystals.p", "rb") as f:
        crytals = pickle.load(f)
    with open(structs_pickle_dir + "/save_reactivities.p", "rb") as f:
        reactivities = pickle.load(f)

if not cap_rhos:
    rho_midpoints = [-1.0]

stationary_args = [reactivities, crystals, sample_n, react_rhos, structs_pickle_dir, output_dir, out_sub_dir, outname, best_dist_func, scaling_func, cap_rhos, shape_slope, shape_intercept]
if '--arg_slice' in opts:
    rm_cv_w = [(float(a) for a in opts["--arg_slice"][1:-1].split(","))]
    print "rm_c_w arg_slice: " + str(rm_cv_w)
else:
    rm_cv_w = zip(OSU.ncycles(rho_midpoints, len(constrain_vals)*len(weights)), cycle(OSU.repeat_elements(constrain_vals, len(rho_midpoints))), cycle(OSU.repeat_elements(weights, len(rho_midpoints)*len(constrain_vals))))

training_results = {}
training_res_dir = OSU.create_directory(output_dir + "/training_results/")
sub_proc_sh_dir = OSU.create_directory(output_dir + "/sub_proc_sh/")

# Option in case processes were aborted
if "--restart" in opts:
    finished = glob.glob(training_res_dir + "save_training_results_*.p")
    for f in finished:
        fname = re.findall("save_training_results_(.+).p$", f)[0]
        f_params = tuple([float(a) for a in fname.split("_")])
        if f_params in rm_cv_w:
            rm_cv_w.remove(f_params)

args_pool = [list(a) + stationary_args for a in rm_cv_w]
remaining_params = len(rm_cv_w)

if num_proc == 1 and not cluster_flag and not load_results:
    # no parallelization on args level
    print "Hit num_proc == 1 and not cluster_flag and not load_results"
    for args in args_pool:
        p_result = PAU.load_train_model_helper(args)[0]
        c = p_result[1]
        if c == []:
            c = ""
        k = CRW_pair(c=c, mr=p_result[0], w=p_result[2])
        training_results[k] = (p_result[3], p_result[4])
        if sub_proc:
            pickle.dump(training_results, open(training_res_dir + "save_training_results_%s_%s_%s.p" % (k.mr, k.c, k.w), "wb"))
elif num_proc == 1 and cluster_flag and not load_results:  # This case is the first executed for the parallel version that utilizes the full cluster.
    # Surrounded job execution code to catch any subproc that doesn't finish to the pickling step.
    # Also acts as a limiter into the number of jobs that can be submitted to the queue at once.
    max_jobs = 511
    jobs_available = 511
    while len(rm_cv_w) > 0:
        sub_proc_dir = OSU.create_directory(output_dir + "sub_proc_out/")
        for param in rm_cv_w:
            param_string = "_".join([str(s) for s in param])
            job_name_param = "_".join([job_name, param_string])[:31]  # Job name can only be up to 31 characters long
            if not OSU.check_file_exists("%snbs_script_%s.sh" % (sub_proc_sh_dir, param_string)):
                header = "##NBS-stdout:%s\n##NBS-stderr:%s\n##NBS-queue:batch\n##NBS-name:\"%s\"\n##NBS-jcoll:\"%s\"\n\nrm %s %s\n" % (sub_proc_dir + job_name_param+".out", sub_proc_dir+job_name_param+".err", job_name_param, job_name, sub_proc_dir + job_name_param+".out", sub_proc_dir+job_name_param+".err")
                OSU.system_command("echo \"%s/usr/bin/time /fs/home/amy35/tools/anaconda/bin/python find_parameters.py -r \'%s\' -c \'%s\' -o %s %s -n %s -p 1 --scaling_func %s --cluster_flag False --sub_proc True --arg_slice \'%s\' --job_name %s --load_results \'False\' --generate_structs \'False\' --cap_rhos %s --structs_pickle_dir %s\"> %snbs_script_%s.sh" % (header, opts['-r'], opts['-c'], opts['-o'], sampling_opts_string, opts['-n'], opts['--scaling_func'], param, job_name_param, cap_rhos, structs_pickle_dir, sub_proc_sh_dir, param_string))
            if jobs_available > 0 and not PAU.check_job_on_queue(job_name_param) and not OSU.check_file_exists("".join([training_res_dir, "save_training_results_", param_string, ".p"])):
                print "/opt/voyager/nbs/bin/jsub %snbs_script_%s.sh -name %s -stdout %snbs_script_%s.out -stderr %snbs_script_%s.err" % (sub_proc_sh_dir, param_string, job_name_param, sub_proc_dir, param_string, sub_proc_dir, param_string)
                OSU.system_command("/opt/voyager/nbs/bin/jsub %snbs_script_%s.sh -name %s -stdout %snbs_script_%s.out -stderr %snbs_script_%s.err" % (sub_proc_sh_dir, param_string, job_name_param, sub_proc_dir, param_string, sub_proc_dir, param_string))
                jobs_available -= 1
            else:
                break
        jobs_available = PAU.wait_jcoll_finish_any(job_name, sub_proc_dir, max_jobs, 60)
        rm_cv_w = [bm_param_set for bm_param_set in rm_cv_w if not OSU.check_file_exists("".join([training_res_dir, "save_training_results_", "_".join([str(b) for b in bm_param_set]), ".p"]))]
        OSU.system_command("echo \"len rm_cv_w: %s\n\" >> %sjcoll_waiting.txt" % (len(rm_cv_w), sub_proc_dir))
        if remaining_params == len(rm_cv_w):
            # send on first in rm_cv_w if there is no change in the number of parameters left
            OSU.system_command("echo \"len rm_cv_w same as last time: %s\n\" >> %sjcoll_waiting.txt" % (rm_cv_w[0], sub_proc_dir))
            param_string = "_".join([str(ele) for ele in rm_cv_w[0]])
            job_name_param = "_".join([job_name, param_string])[:31]  # job name can only be 1-31 characters long
            OSU.system_command('echo "/opt/voyager/nbs/bin/jsub %snbs_script_%s.sh -name %s -stdout %snbs_script_%s.out -stderr %snbs_script_%s.err" >> %sjcoll_waiting.txt' % (sub_proc_sh_dir, param_string, job_name_param, sub_proc_dir, param_string, sub_proc_dir, param_string, sub_proc_dir))
            OSU.system_command("/opt/voyager/nbs/bin/jsub %snbs_script_%s.sh -name %s -stdout %snbs_script_%s.out -stderr %snbs_script_%s.err" % (sub_proc_sh_dir, param_string, job_name_param, sub_proc_dir, param_string, sub_proc_dir, param_string))
            OSU.system_command("echo \"waiting for this output file: %s\n%s\n\" >> %sjcoll_waiting.txt" % (rm_cv_w[0], "".join([training_res_dir, "save_training_results_", "_".join([str(b) for b in rm_cv_w[0]]), ".p"]), sub_proc_dir))
            PAU.wait_for_jcoll(job_name, sub_proc_dir, 30)
            OSU.system_command("echo \"len rm_cv_w same as last time: %s\n\" >> %sjcoll_waiting.txt" % ("".join([training_res_dir, "save_training_results_", "_".join([str(b) for b in rm_cv_w[0]]), ".p"]), sub_proc_dir))
            if OSU.check_file_exists("".join([training_res_dir, "save_training_results_", "_".join([str(b) for b in rm_cv_w[0]]), ".p"])):
                rm_cv_w.pop(0)
                OSU.system_command("echo \"len rm_cv_w after running one: %s\n\" >> %sjcoll_waiting.txt" % (len(rm_cv_w), sub_proc_dir))
            else:
                OSU.system_exit("Did not find %s" % ("".join([training_res_dir, "save_training_results_", "_".join([str(b) for b in rm_cv_w[0]]), ".p"])))

        remaining_params = len(rm_cv_w)
        print "remaining_params: " + str(remaining_params)

    training_res_pickles = OSU.load_all_pickles(training_res_dir)
    for pickle_value in training_res_pickles.values():
        training_results.update(pickle_value)
    max_F = max([f for dr, f in training_results.values()])
    max_F_keys = [k for k in sorted(training_results) if training_results[k][1] == max_F]
    with open(output_dir + "/best_params.txt", "w") as f:
        f.write(str(max_F) + " : Sum of F-scores\n")
        f.write(str(max_F_keys) + " : Constraint values of best sum F-Scores\n")
else:
    raise Exception("Cluster case not implemented yet")

if load_results and not sub_proc:
    if num_proc != 1:
        training_res_pickles = OSU.load_all_pickles(training_res_dir)
        for pickle_value in training_res_pickles.values():
            training_results.update(pickle_value)
    print "training_results.values(): " + str(training_results.values())
    max_F = max([f for dr, f in training_results.values()])
    max_F_keys = [k for k in sorted(training_results) if training_results[k][1] == max_F]
    with open(output_dir + "/best_params.txt", "w") as f:
        f.write(str(max_F) + " : Sum of F-scores\n")
        f.write(str(max_F_keys) + " : Constraint values of best sum F-Scores\n")
