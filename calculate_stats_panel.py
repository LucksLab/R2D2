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

# setup environment variables
LucksLabUtils_config.config("Quest_R2D2")
opts = OSU.getopts("o:c:r:p:", ["shape_intercept=", "shape_slope="])
print opts

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

for k, v in reactivities.iteritems():
    react_rhos = v[0]
    react_seq = v[1]
    ck = k.split('-')[0]  # corresponding crystal key from the reactivity key
    cryst_seq = crystals[ck][1]
    ind_of_match = react_seq.index(cryst_seq)
    end_of_match = ind_of_match + len(cryst_seq)
    # renormalize rhos based on matching regions
    react_rhos = SU.recalc_rhos(react_rhos, ind_of_match, end_of_match)

    # SHAPE folding
    reactivities[k][0] = react_rhos
    reactivities[k][1] = react_seq[ind_of_match:end_of_match]
    SU.rhos_list_to_file(react_rhos, reactivities[k][2]+".rho")
    seqfile = NAU.make_seq(reactivities[k][1], reactivities[k][2]+".seq")
    SU.runRNAstructure_fold(seqfile, reactivities[k][2]+".ct", shapefile=reactivities[k][2]+".rho", p=num_proc, shape_intercept=shape_intercept, shape_slope=shape_slope)
    SU.runRNAstructure_CircleCompare(reactivities[k][2]+".ct", crystals[ck][3], reactivities[k][2]+".ps")
    OSU.system_command("convert %s.ps %s.jpg" % (reactivities[k][2], reactivities[k][2]))
    with open(reactivities[k][2]+".stats", "w") as f:
        react_mat = SU.ct_struct_to_binary_mat(SU.get_ct_structs(reactivities[k][2]+".ct"))[0]
        f.write(str(SU.calc_benchmark_statistics_matrix(react_mat, crystals[ck][2])))

    # no SHAPE folding
    SU.runRNAstructure_fold(seqfile, reactivities[k][2]+"_noshape.ct", p=num_proc, shape_intercept=shape_intercept, shape_slope=shape_slope)
    SU.runRNAstructure_CircleCompare(reactivities[k][2]+"_noshape.ct", crystals[ck][3], reactivities[k][2]+"_noshape.ps")
    OSU.system_command("convert %s_noshape.ps %s_noshape.jpg" % (reactivities[k][2], reactivities[k][2]))
    with open(reactivities[k][2]+"_noshape.stats", "w") as f:
        react_mat = SU.ct_struct_to_binary_mat(SU.get_ct_structs(reactivities[k][2]+"_noshape.ct"))[0]
        f.write(str(SU.calc_benchmark_statistics_matrix(react_mat, crystals[ck][2])))

    # R2D2
    # taking code from cotranscriptional case (PCSU.run_cotrans_length) which has extraneous parts in this use case
    # few lines in PCSU.run_cotrans_length made it unable to be used for this case. ex. length_key
    e = 50000
    #output_file_prefix = re.findall("([^/]+).rnm$", rnmfile)[0]
    scaling_fns = {"D": SU.invert_scale_rho_vec, "U": SU.scale_vec_avg1, "K": SU.cap_rho_or_ct_list}
    scaling_func = "K"
    scale_rho_max = 1.0
    constrained_c = 3.5
    cap_rhos = True
    weight_paired = 0.8
    draw_dir = OSU.create_directory("%s/%s" % (R2D2_output_dir, k))
    bm_R2D2_all = defaultdict(list)

    for rnum in range(1,101):
        sampled_structs_count = defaultdict(int)
        sampled_structs = set()

        # Vanilla Sampling
        structs, structs_labels = SU.RNAstructure_sample(reactivities[k][2], e, R2D2_output_dir, label="noshape", num_proc=1, wn_tag=k)
        sampled_structs.update(structs_labels)
        OSU.increment_dict_counts(sampled_structs_count, structs)

        # Sampling with SHAPE constraints
        structs, structs_labels = SU.RNAstructure_sample(reactivities[k][2], e, R2D2_output_dir, shapefile=reactivities[k][2]+".rho", label="shape", num_proc=1, wn_tag=k)
        sampled_structs.update(structs_labels)
        OSU.increment_dict_counts(sampled_structs_count, structs)

        # Sampling with hard constraints
        XB = SU.get_indices_rho_gt_c(react_rhos, constrained_c, one_index=True)  # RNAstructure is 1-indexed
        SU.make_constraint_file(reactivities[k][2]+".con", [], XB, [], [], [], [], [], [])
        structs, structs_labels = SU.RNAstructure_sample(reactivities[k][2], e, R2D2_output_dir, constraintfile=reactivities[k][2]+".con", label="constrained_"+str(constrained_c), num_proc=1, wn_tag=k)
        sampled_structs.update(structs_labels)
        OSU.increment_dict_counts(sampled_structs_count, structs)

        # Compressing sampled structures further by removing duplicates sampled by multiple methods. Keeping track of this though.
        # Saving more than I need to in this use case... ex. energies
        sampled_structs = SU.merge_labels(sampled_structs, to_string=False)
        structs = [t[0].split(",") for t in sampled_structs]
        SU.cts_to_file(structs, cryst_seq, reactivities[k][2]+"_unique.ct")
        SU.runRNAstructure_efn2(reactivities[k][2]+"_unique.ct", reactivities[k][2] + ".efn2")
        free_energies = SU.get_free_energy_efn2(reactivities[k][2] + ".efn2")

        if cap_rhos:
            scaled_rhos = scaling_fns[scaling_func](react_rhos, scale_rho_max)
        else:
            scaled_rhos = scaling_fns[scaling_func](react_rhos)
        with open(reactivities[k][2]+".best_scaled_rho", "w") as f:
            f.write("\n".join(["\t".join([str(zi), str(zr)]) for zi,zr in enumerate(scaled_rhos)]))

        # Compute distances between scaled rhos and paired-vectors from drawn structures
        binary_structs = SU.ct_struct_to_binary_vec(structs)
        distances = []
        for s in binary_structs:
            distances.append(SU.calc_bp_distance_vector_weighted(s, scaled_rhos, scaling_func=scaling_func, invert_struct="D" != scaling_func, paired_weight=weight_paired))
        min_distance = min(distances)
        min_dist_indices = [i for i, v in enumerate(distances) if v == min_distance]

        # compare R2D2 against crystal structure
        for mdi in min_dist_indices:
            react_mat = SU.ct_struct_to_binary_mat(structs[mdi])
            curr_prefix = "%s_%s_R2D2"%(reactivities[k][2], mdi)
            curr_stats = SU.calc_benchmark_statistics_matrix(react_mat, crystals[ck][2])
            with open(curr_prefix + ".stats", "w") as f:
                f.write(str(curr_stats))
            for csk in curr_stats:
                bm_R2D2_all[csk].append(curr_stats[csk])
            #make file
            SU.ct_list_to_file(structs[mdi], cryst_seq, curr_prefix+".ct")
            SU.runRNAstructure_CircleCompare(curr_prefix+".ct", crystals[ck][3], curr_prefix+".ps")
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
            SU.run_ct2dot(reactivities[k][2]+"_unique.ct", seen_snum[sf], draw_outname_pre + ".dbn")
    # finished R2D2 iterations for this dataset
    # report all and avg R2D2 results
    with open("%s_R2D2_all.stats" % (reactivities[k][2]), "w") as f:
        f.write(str(bm_R2D2_all) + "\n")
        for key in bm_R2D2_all:
            f.write("%s\t%s\n" % (key, sum(bm_R2D2_all[key]) / len(bm_R2D2_all[key])))
