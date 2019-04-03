"""
Calculate the performance statistics of RNAstructure-Fold on individual replicates in benchmark panel.

Options:
-r                = reactivities files regex
-c                = crystal structure (.ct) files regex
-o                = output directory
-p                = number of threads to use, default 1
--shape_intercept = Intercept used with SHAPE restraints in RNAstructure, default -0.3
--shape_slope     = Slope used with SHAPE restraints in RNAstructure, default 1.1

Version: 0.0.1
Author: Angela M Yu, 2014-2019

Copyright (C) 2014-2019  Julius B. Lucks and Angela M Yu.
All rights reserved.
Distributed under the terms of the GNU General Public License, see 'LICENSE'.
"""

from __future__ import division  # allows division of ints to output decimals
import SU
import OSU
import PAU
import NAU
import LucksLabUtils_config
import glob
import cPickle as pickle
from collections import defaultdict
import re
from itertools import repeat
from multiprocessing import Pool, Lock
import numpy
from sys import maxsize
from sklearn.decomposition import PCA
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt


# setup environment variables
LucksLabUtils_config.config("Quest_R2D2")
opts = OSU.getopts("o:c:r:p:", ["shape_intercept=", "shape_slope="])
print opts
numpy.set_printoptions(threshold=maxsize)
plt.style.use('seaborn-whitegrid')

reactivities_files = glob.glob(opts['-r'])
crystal_files = glob.glob(opts['-c'])
output_dir = OSU.create_directory(opts['-o'])
R2D2_output_dir = OSU.create_directory(output_dir + "/R2D2_intermediate/")
num_proc = int(opts['-p']) if '-p' in opts else 1
shape_intercept = float(opts["--shape_intercept"]) \
                    if "--shape_intercept" in opts else -0.3
shape_slope = float(opts["--shape_slope"]) if "--shape_slope" in opts else 1.1

reactivities = PAU.parse_input_panels(reactivities_files, output_dir)
crystals = PAU.parse_input_panels(crystal_files, output_dir)

if set(crystals.keys()) != set([rk.split('-')[0] for rk in reactivities.keys()]):
    raise Exception("Keys of crystal structures and reactivities not equal")

pickle.dump(crystals, open(output_dir + "/crystals.p", "wb"))

print crystals.keys()
print reactivities.keys()

def R2D2_process_wrapper(arg_tuple):
    return R2D2_process(*arg_tuple)


def R2D2_process(input_prefix, R2D2_output_dir, draw_dir, react_rhos, crystals_ck, rnum, parallel=False):
    """
    Slightly reduced version of a R2D2 process for this benchmarking code.
    # taking code from cotranscriptional case (PCSU.run_cotrans_length) which has extraneous parts in this use case
    # few lines in PCSU.run_cotrans_length made it unable to be used for this case. ex. length_key

    Options:
    input_prefix - full path plus input prefix of reactivities
    R2D2_output_dir - R2D2 output directoory
    draw_dir - directory for .dbn output
    react_rhos - rho reactivities
    rnum - iteration number, names output files accordingly
    """
    e = 50000
    fname = re.findall("([^/]+)$", input_prefix)[0]
    output_prefix = "%s/%s_%s" % (R2D2_output_dir, fname, rnum)
    scaling_fns = {"D": SU.invert_scale_rho_vec, "U": SU.scale_vec_avg1, "K": SU.cap_rho_or_ct_list}
    scaling_func = "K"
    scale_rho_max = 1.0
    constrained_c = 3.5
    cap_rhos = True
    weight_paired = 0.8

    sampled_structs_count = defaultdict(int)
    sampled_structs = set()

    # Vanilla Sampling
    structs, structs_labels = SU.RNAstructure_sample(input_prefix, e, R2D2_output_dir, label="noshape", num_proc=1, wn_tag="_%s"%(rnum), lock=lock)
    sampled_structs.update(structs_labels)
    OSU.increment_dict_counts(sampled_structs_count, structs)

    # Sampling with SHAPE constraints
    structs, structs_labels = SU.RNAstructure_sample(input_prefix, e, R2D2_output_dir, shapefile=input_prefix+".rho", label="shape", num_proc=1, wn_tag="_%s"%(rnum), lock=lock)
    sampled_structs.update(structs_labels)
    OSU.increment_dict_counts(sampled_structs_count, structs)

    # Sampling with hard constraints
    XB = SU.get_indices_rho_gt_c(react_rhos, constrained_c, one_index=True)  # RNAstructure is 1-indexed
    SU.make_constraint_file(output_prefix+".con", [], XB, [], [], [], [], [], [])
    structs, structs_labels = SU.RNAstructure_sample(input_prefix, e, R2D2_output_dir, constraintfile=output_prefix+".con", label="constrained_"+str(constrained_c), num_proc=1, wn_tag="_%s"%(rnum), lock=lock)
    sampled_structs.update(structs_labels)
    OSU.increment_dict_counts(sampled_structs_count, structs)

    # Compressing sampled structures further by removing duplicates sampled by multiple methods. Keeping track of this though.
    # Saving more than I need to in this use case... ex. energies
    sampled_structs = SU.merge_labels(sampled_structs, to_string=False)
    structs = [t[0].split(",") for t in sampled_structs]
    SU.cts_to_file(structs, crystals_ck[1], output_prefix+"_unique.ct")
    SU.runRNAstructure_efn2(output_prefix+"_unique.ct", output_prefix + ".efn2")
    free_energies = SU.get_free_energy_efn2(output_prefix + ".efn2")

    if cap_rhos:
        scaled_rhos = scaling_fns[scaling_func](react_rhos, scale_rho_max)
    else:
        scaled_rhos = scaling_fns[scaling_func](react_rhos)
    with open(input_prefix+".best_scaled_rho", "w") as f:
        f.write("\n".join(["\t".join([str(zi), str(zr)]) for zi,zr in enumerate(scaled_rhos)]))

    # Compute distances between scaled rhos and paired-vectors from drawn structures
    binary_structs = SU.ct_struct_to_binary_vec(structs)
    distances = []
    for s in binary_structs:
        distances.append(SU.calc_bp_distance_vector_weighted(s, scaled_rhos, scaling_func=scaling_func, invert_struct="D" != scaling_func, paired_weight=weight_paired))
    min_distance = min(distances)
    min_dist_indices = [i for i, v in enumerate(distances) if v == min_distance]

    # compare R2D2 against crystal structure
    selected_react_mats = []
    for mdi in min_dist_indices:
        react_mat = SU.ct_struct_to_binary_mat(structs[mdi])
        selected_react_mats.append(numpy.matrix(react_mat))
        curr_prefix = "%s_%s_R2D2"%(output_prefix, mdi)
        curr_stats = SU.calc_benchmark_statistics_matrix(react_mat, crystals_ck[2])
        with open(curr_prefix + ".stats", "w") as f:
            f.write(str(curr_stats))
        #make file
        SU.ct_list_to_file(structs[mdi], crystals_ck[1], curr_prefix+".ct")
        SU.runRNAstructure_CircleCompare(curr_prefix+".ct", crystals_ck[3], curr_prefix+".ps")
        OSU.system_command("convert %s.ps %s.jpg" % (curr_prefix, curr_prefix))

    # saving R2D2 results
    R2D2_save = {}
    R2D2_save["structs"] = structs
    R2D2_save["distances"] = distances
    R2D2_save["min_dist_indices"] = min_dist_indices
    R2D2_save["min_distance"] = min_distance
    R2D2_save["scaled_rhos"] = scaled_rhos
    R2D2_save["react_mat"] = react_mat
    pickle.dump(R2D2_save, open(curr_prefix + ".p", "wb"))

    # output .dbn's like in normal R2D2 process
    # code taken from PCSU.generate_best_struct_images
    # PCSU.generate_best_struct_images contained some extraneuous calls for this use case. ex. draw_all = False, running VARNA
    seen_snum = []
    iter_dbn_dir = OSU.create_directory("%s/%s" % (draw_dir, rnum))
    for snum in min_dist_indices:
        seen_snum.append(snum)
    for sf in range(len(seen_snum)):
        draw_outname_pre = "%s/%snt_%s" % (iter_dbn_dir, len(react_rhos), seen_snum[sf])
        if len(seen_snum) > 1:
            draw_outname_pre += "_mult" + str(sf)
        SU.run_ct2dot(output_prefix+"_unique.ct", seen_snum[sf], draw_outname_pre + ".dbn")

    # return curr_stats and selected structures
    return curr_stats, selected_react_mats


def parse_R2D2_process_output(out_stats, selected_react_mats, bm_R2D2_all, R2D2_all_selected, num_selected, R2D2_selected_binary_cts):
    """ parse R2D2 process output """
    for csk, value in out_stats.items():
        bm_R2D2_all[csk].append(value)
    for srm in selected_react_mats:
        R2D2_all_selected = numpy.add(R2D2_all_selected, srm)
        num_selected += 1
        R2D2_selected_binary_cts.append(SU.binary_mat_to_binary_ct(srm))
    return bm_R2D2_all, R2D2_all_selected, num_selected, R2D2_selected_binary_cts


def init(l):
    """
    RNAstructure exe's look like they collide when run in parallel, so making a lock
    """
    global lock
    lock = l


def write_reactivities_in_ct(ct, react, outfile):
    react_ct_zip = zip(react, [int(c) for c in ct])
    with open(outfile, "w") as f:
        f.write(str(react_ct_zip) + "\n")
        f.write("Paired:\t" + "\t".join([str(r) for r, c in react_ct_zip if c > 0]) + "\n")
        f.write("Unpaired:\t" + "\t".join([str(r) for r, c in react_ct_zip if c == 0]) + "\n")
        f.write("Avg paired:\t" + str(sum([r for r, c in react_ct_zip if c > 0]) / len(react_ct_zip)) + "\n")
        f.write("Avg unpaired:\t" + str(sum([r for r, c in react_ct_zip if c == 0]) / len(react_ct_zip)) + "\n")


if num_proc > 1:
    l = Lock()
    pool = Pool(processes=num_proc, initializer=init, initargs=(l,))

for k, v in reactivities.iteritems():
    react_rhos = v[0]
    react_seq = v[1]
    ck = k.split('-')[0]  # corresponding crystal key from the reactivity key
    cryst_seq = crystals[ck][1]
    ind_of_match = react_seq.index(cryst_seq)
    end_of_match = ind_of_match + len(cryst_seq)
    # renormalize rhos based on matching regions
    react_rhos = SU.recalc_rhos(react_rhos, ind_of_match, end_of_match)
    reactivities[k][0] = react_rhos
    reactivities[k][1] = react_seq[ind_of_match:end_of_match]

    # SHAPE folding
    SU.rhos_list_to_file(react_rhos, reactivities[k][2]+".rho")
    seqfile = NAU.make_seq(reactivities[k][1], reactivities[k][2]+".seq")
    SU.runRNAstructure_fold(seqfile, reactivities[k][2]+"_Fold_shape.ct", shapefile=reactivities[k][2]+".rho", p=num_proc, shape_intercept=shape_intercept, shape_slope=shape_slope)
    SU.runRNAstructure_CircleCompare(reactivities[k][2]+"_Fold_shape.ct", crystals[ck][3], reactivities[k][2]+"_Fold_shape.ps")
    OSU.system_command("convert %s_Fold_shape.ps %s_Fold_shape.jpg" % (reactivities[k][2], reactivities[k][2]))
    with open(reactivities[k][2]+".stats", "w") as f:
        fold_shape_ct = SU.get_ct_structs(reactivities[k][2]+"_Fold_shape.ct")[0]
        react_mat = SU.ct_struct_to_binary_mat(fold_shape_ct)
        f.write(str(SU.calc_benchmark_statistics_matrix(react_mat, crystals[ck][2])))
    write_reactivities_in_ct(fold_shape_ct, react_rhos, reactivities[k][2]+"_Fold_shape_ct_react.txt")

    # no SHAPE folding
    SU.runRNAstructure_fold(seqfile, reactivities[k][2]+"_Fold_noshape.ct", p=num_proc, shape_intercept=shape_intercept, shape_slope=shape_slope)
    SU.runRNAstructure_CircleCompare(reactivities[k][2]+"_Fold_noshape.ct", crystals[ck][3], reactivities[k][2]+"_Fold_noshape.ps")
    OSU.system_command("convert %s_Fold_noshape.ps %s_Fold_noshape.jpg" % (reactivities[k][2], reactivities[k][2]))
    with open(reactivities[k][2]+"_Fold_noshape.stats", "w") as f:
        fold_noshape_ct = SU.get_ct_structs(reactivities[k][2]+"_Fold_noshape.ct")[0]
        react_mat = SU.ct_struct_to_binary_mat(fold_noshape_ct)
        f.write(str(SU.calc_benchmark_statistics_matrix(react_mat, crystals[ck][2])))
    write_reactivities_in_ct(fold_noshape_ct, react_rhos, reactivities[k][2]+"_Fold_noshape_ct_react.txt")

    # R2D2
    # taking code from cotranscriptional case (PCSU.run_cotrans_length) which has extraneous parts in this use case
    # few lines in PCSU.run_cotrans_length made it unable to be used for this case. ex. length_key
    draw_dir = OSU.create_directory("%s/%s" % (R2D2_output_dir, k))
    bm_R2D2_all = defaultdict(list)
    R2D2_all_selected = numpy.zeros((len(cryst_seq), len(cryst_seq)))
    R2D2_selected_binary_cts = []
    num_selected = 0
    if num_proc > 1:
        args_pool = zip(repeat(reactivities[k][2]), repeat(R2D2_output_dir), repeat(draw_dir), repeat(react_rhos), repeat(crystals[ck]), range(100), repeat(True))
        for out_stats, selected_react_mats in pool.imap(R2D2_process_wrapper, args_pool):
            bm_R2D2_all, R2D2_all_selected, num_selected, R2D2_selected_binary_cts = parse_R2D2_process_output(out_stats, selected_react_mats, bm_R2D2_all, R2D2_all_selected, num_selected, R2D2_selected_binary_cts)
        del args_pool
    else:
        for rnum in range(1,101):
            out_stats, selected_react_mats = R2D2_process(input_prefix=reactivities[k][2], R2D2_output_dir=R2D2_output_dir, draw_dir=draw_dir, rnum=rnum)
            bm_R2D2_all, R2D2_all_selected, num_selected, R2D2_selected_binary_cts = parse_R2D2_process_output(out_stats, selected_react_mats, bm_R2D2_all, R2D2_all_selected, num_selected, R2D2_selected_binary_cts)

    # finished R2D2 iterations for this dataset
    # report all and avg R2D2 results
    with open("%s_R2D2_all.stats" % (reactivities[k][2]), "w") as f:
        f.write(str(bm_R2D2_all) + "\n")
        for key in bm_R2D2_all:
            f.write("%s\t%s\n" % (key, sum(bm_R2D2_all[key]) / len(bm_R2D2_all[key])))

    # find basepairs that occurred over 50% of the time, then calculate accuracy
    R2D2_pairs = R2D2_all_selected / num_selected
    with open("%s_R2D2_pairs.txt" % (reactivities[k][2]), "w") as f:
        f.write("\n".join(["\t".join([str(bp) for bp in row]) for row in R2D2_pairs.tolist()]) + "\n")
    R2D2_consensus = R2D2_pairs.round()
    with open("%s_R2D2_consensus.txt" % (reactivities[k][2]), "w") as f:
        f.write("\n".join(["\t".join([str(bp) for bp in row]) for row in R2D2_consensus.tolist()]) + "\n")
    with open("%s_R2D2_consensus.stats" % (reactivities[k][2]), "w") as f:
        f.write(str(SU.calc_benchmark_statistics_matrix(R2D2_consensus, crystals[ck][2])))
    write_reactivities_in_ct(SU.binary_mat_to_binary_ct(R2D2_consensus), react_rhos, reactivities[k][2]+"_R2D2_consensus_ct_react.txt")
    R2D2_consensus_ct = SU.binary_mat_to_ct(R2D2_consensus)
    SU.ct_list_to_file(R2D2_consensus_ct, cryst_seq, "%s_R2D2_consensus.ct" % (reactivities[k][2]))

    # reactivities compared to crystal structure
    write_reactivities_in_ct(crystals[ck][4], react_rhos, reactivities[k][2]+"_crystal_ct_react.txt")

    # PCA
    pca = PCA(n_components=2)
    X = SU.ct_struct_to_binary_vec([fold_shape_ct, fold_noshape_ct, R2D2_consensus_ct, crystals[ck][4]])
    X += R2D2_selected_binary_cts
    Y = ["Fold_SHAPE", "Fold_noSHAPE", "R2D2_consensus", "crystal"]
    Y += ["R2D2_iteration"] * len(R2D2_selected_binary_cts)
    colors = ['r', 'g', 'b', 'm'] + ['c'] * len(R2D2_selected_binary_cts)
    principal_components = pca.fit_transform(numpy.array(X))

    with open("%s_PCA_coords.txt" % (reactivities[k][2]), "w") as f:
        f.write("\n".join(["\t".join([str(coord) for coord in row]) for row in principal_components.tolist()]) + "\n")

    fig = plt.figure(figsize = (8,8))
    ax = fig.add_subplot(1,1,1) 
    ax.set_xlabel('Principal Component 1', fontsize = 15)
    ax.set_ylabel('Principal Component 2', fontsize = 15)
    ax.set_title('PCA', fontsize = 20)
    ax.scatter(principal_components[:,0], principal_components[:,1], c=colors, s=23)
    # adding individual points separately to get legend to work, and maintain correct axes window
    for i in range(5):
        ax.scatter(principal_components[i,0], principal_components[i,1], c=colors[i], s=23, label=Y[i])
    ax.legend()
    plt.savefig("%s_PCA.png" % (reactivities[k][2]))

    del X, principal_components
