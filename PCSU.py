"""
Parallel Cotranscriptional Structure Utilities (PCSU)

Version: 0.0.1
Author: Angela M Yu, 2014-2016

Copyright (C) 2016  Julius B. Lucks, Angela M Yu.
All rights reserved.
Distributed under the terms of the GNU General Public License, see 'LICENSE'.
"""

import SU
import OSU
import re
from collections import defaultdict
from itertools import cycle
import cPickle as pickle
import VIU
from multiprocessing import Pool


def calculate_function_helper(args):
    """
    Helper function used to run the function func with arguements args
    """
    return calculate_function(*args)


def calculate_function(func, args):
    """
    Used to run the function func with arguements args
    """
    return func(*args)


def run_cotrans_length_helper(args):
    """
    Helper function used to run the analysis of individual lengths
    """
    return run_cotrans_length(*args)


def run_output_multiprocessing_pool(function, args_pool, num_proc):
    """
    Control multiproccessing of a function with many sets of arguments as defined in args_pool.
    num_proc: number of proccesses to use
    """
    for i in xrange(0, len(args_pool), num_proc):
        print "running output args_slice nums: " + str(i)
        args_slice = args_pool[i:(i+num_proc)]
        pool = Pool(processes=num_proc)
        pool.imap_unordered(function, args_slice)
        pool.close()
        pool.join()
    return


def run_cotrans_length(file_l, output_dir, ct_dir, pickle_dir, adapterseq, endcut, pol_fp, e, constrained_c, cap_rhos, scale_rho_max, scaling_func, weight_paired):
    """
    Main analysis of a length in the CoTranscriptional SHAPE-Seq experiment.

    file_l: 	reactivity file at one length
    output_dir: output directory
    ct_dir: 	directory to output .ct files
    pickle_dir: directory to output pickles
    adapterseq: adapter sequence to be trimmed
    endcut: 	Number of bases to be ommitted from calculations from the 3' end of the length.
                Ex. Ignoring the last base would correspond to endcut=-1
    pol_fp:     Number of bases occupied by polymerase. Globally sets for each length.
                    Ex. Polymerase foot print of size 15 => pol_fp=15
    e: 		Number of structures sampled by each sampling method
    constrained_c: Any rho value greater than or equal to this value will be forced as unpaired when sampling with hard constraints
    cap_rhos: 	Flag to have a max cutoff when calculating distances for choosing the best structure
    scale_rho_max: If the cap_rhos flag is True or the 'D' distance is used, this value is the max value and all values greater than it are set to this max value
    scaling_func: Choice of distance function when choosing the best structure:
                D: Make reactivities to be bound between [0,1]
                U: Rescale structures to average to 1
                K: Scaling of reactivities and structures are kept
    weight_paired: Weight given to paired regions in distance calculations
    """

    # Input handling
    scaling_fns = {"D": SU.invert_scale_rho_vec, "U": SU.scale_vec_avg1, "K": SU.cap_rho_or_ct_list}
    fname = re.findall("([^/]+).txt$", file_l)
    output_file_prefix = output_dir+"/"+fname[0]

    # theta, rho and seq have been cut by -adapter_len + endcut - pol_fp
    pos, rho_full, theta, rho, seq, rc_flag, rc_sum, rc_untreated_sum, rc_treated_sum = SU.parse_reactivity_rho(file_l, adapterseq, output_file_prefix, endcut - pol_fp)

    length_key = len(pos) + abs(endcut - pol_fp)
    file_data_length_key = {}
    file_data_length_key["filename"] = fname[0]
    file_data_length_key["rho"] = rho
    file_data_length_key["seq"] = seq
    del theta, rc_treated_sum, rc_untreated_sum, rc_sum

    sampled_structs_count = defaultdict(int)
    sampled_structs = set()

    # Vanilla Sampling
    structs, structs_labels = SU.RNAstructure_sample(output_file_prefix, e, output_dir, label="noshape", num_proc=1, wn_tag=str(length_key))
    sampled_structs.update(structs_labels)
    OSU.increment_dict_counts(sampled_structs_count, structs)

    # Sampling with SHAPE constraints
    structs, structs_labels = SU.RNAstructure_sample(output_file_prefix, e, output_dir, shapefile=output_file_prefix+".rho", label="shape", num_proc=1, wn_tag=str(length_key))
    sampled_structs.update(structs_labels)
    OSU.increment_dict_counts(sampled_structs_count, structs)

    # Sampling with hard constraints
    XB = SU.get_indices_rho_gt_c(rho, constrained_c, one_index=True)  # RNAstructure is 1-indexed
    SU.make_constraint_file(output_file_prefix+".con", [], XB, [], [], [], [], [], [])
    structs, structs_labels = SU.RNAstructure_sample(output_file_prefix, e, output_dir, constraintfile=output_file_prefix+".con", label="constrained_"+str(constrained_c), num_proc=1, wn_tag=str(length_key))
    sampled_structs.update(structs_labels)
    OSU.increment_dict_counts(sampled_structs_count, structs)

    # Compressing sampled structures further by removing duplicates sampled by multiple methods. Keeping track of this though.
    sampled_structs = SU.merge_labels(sampled_structs, to_string=False)
    structs = [t[0].split(",") for t in sampled_structs]
    SU.cts_to_file(structs, seq, ct_dir+fname[0]+"_unique.ct")
    #JBL- entering debugging here - breakpoint 2 - have checked reactivity parsing, endcutting and renormalization, structure sampling by all three methods
    import ipdb; ipdb.set_trace() 
    SU.runRNAstructure_efn2(ct_dir+fname[0]+"_unique.ct", output_file_prefix + ".efn2")
    free_energies = SU.get_free_energy_efn2(output_file_prefix + ".efn2")
    file_data_length_key["free_energies"] = free_energies
    del free_energies
    file_data_length_key["structs_count"] = sampled_structs_count
    file_data_length_key["structs"] = structs
    file_data_length_key["sampled_structs"] = sampled_structs

    if cap_rhos:
        scaled_rhos = scaling_fns[scaling_func](rho, scale_rho_max)
    else:
        scaled_rhos = scaling_fns[scaling_func](rho)

    file_data_length_key["scaled_rhos"] = scaled_rhos
    with open(output_file_prefix+".best_scaled_rho", "w") as f:
        f.write("\n".join(["\t".join([str(zi), str(zr)]) for zi,zr in enumerate(scaled_rhos)]))

    # Compute distances between scaled rhos and paired-vectors from drawn structures
    binary_structs = SU.ct_struct_to_binary_vec(file_data_length_key["structs"])
    distances = []
    for s in binary_structs:
        distances.append(SU.calc_bp_distance_vector_weighted(s, scaled_rhos, scaling_func=scaling_func, invert_struct="D" != scaling_func, paired_weight=weight_paired))
    file_data_length_key["distances"] = distances
    min_distance = min(distances)
    file_data_length_key["min_dist_indices"] = [i for i, v in enumerate(distances) if v == min_distance]
    file_data_length_key["min_distance"] = min_distance
    file_data_length_key["rc_flag"] = rc_flag
    struct_distances = zip(file_data_length_key["structs"], distances)

    pickle.dump(file_data_length_key, open(pickle_dir + "file_data_" + str(length_key) + ".p", "wb"))
    return length_key, file_data_length_key, struct_distances, len(file_data_length_key["min_dist_indices"]), rho_full, rho


def generate_DG_output(cotrans, start=-1, end=-1):
    """
    Generate DG state plots and .dump file
    """
    if start == -1:
        start = sorted(cotrans.file_data)[0]
    if end == -1:
        end = sorted(cotrans.file_data)[-1]
    print "generate_DG_output: " + str(start) + " " + str(end)
    with open(cotrans.output_dir + "/DG_state_plot.dump", 'w') as dump:
        dump.write("nt\tDG\tmfe_flag\tbest_flag\tdistance\trc_flag\n")
        for length in sorted(cotrans.file_data):
            DG = cotrans.file_data[length]["free_energies"]
            min_DG = min(DG)
            best = cotrans.file_data[length]["min_dist_indices"]  # list of struct_num of min_distance

            line = ["\t".join([str(le[0]),  # nt
                               str(le[1]),  # DG
                               str(int(min_DG == le[1])),  # mfe_flag
                               str(int(c in best)),  # best_flag
                               str(cotrans.file_data[length]["distances"][c]),  # distance
                               str(cotrans.file_data[length]["rc_flag"])])  # rc_flag
                    for le, c in zip(zip(cycle([length]), DG), range(len(DG)))]
            dump.write("\n".join(line))
            dump.write("\n")

    print "R < make_DG_state_plot.R --no-save --args %s/DG_state_plot.pdf %s/DG_state_plot.dump %s %s" % (cotrans.output_dir, cotrans.output_dir, start, end)
    OSU.system_command("R < make_DG_state_plot.R --no-save --args %s/DG_state_plot.pdf %s/DG_state_plot.dump %s %s" % (cotrans.output_dir, cotrans.output_dir, start, end))
    return


def generate_best_struct_images(cotrans_length, length, longest_length, varna_num, zero_padding, draw_dir, ct_dir, draw_all, most_count_tie_break):
    """
    Generates images of the minimum distance structures at a given length
    """
    print "generate_best_struct_images: " + str(length)
    fname = cotrans_length["filename"]
    rho_varna = "\"" + ";".join([str(r) for r in cotrans_length["rho"]]) + "\""

    seen_snum = []
    mult_images = []
    if cotrans_length["rc_flag"]:  # plot only structures with a good RC
        for snum in cotrans_length["min_dist_indices"]:
            seen_snum.append(snum)
    for sf in range(len(seen_snum)):
        draw_outname_pre = "%s/%s_%s_%s" % (draw_dir, fname, seen_snum[sf], str(varna_num).zfill(zero_padding))
        if len(seen_snum) > 1 and draw_all:
            draw_outname_pre += "_mult" + str(sf)
            mult_images.append(draw_outname_pre + "_structure.png")
        elif len(seen_snum) > 1 and not draw_all:
            # draw only the structure with the most supporting counts from the Boltzman samples
            structs_str = [",".join(s) for s in cotrans_length["structs"]]
            if most_count_tie_break:
                tie_break = max(dict((k, v) for k, v in cotrans_length["sampled_structs_count"].iteritems() if structs_str.index(k) in cotrans_length["min_dist_indices"]))
                tie_break_i = structs_str.index(tie_break)
            else:
                tie_break = min([cotrans_length["free_energies"][sn] for sn in cotrans_length["min_dist_indices"]])
                tie_break_i = [cotrans_length["free_energies"][sn] for sn in cotrans_length["min_dist_indices"]].index(tie_break)
                tie_break_i = cotrans_length["min_dist_indices"][tie_break_i]
            if tie_break_i != seen_snum[sf]:
                continue
        rho_varna = rho_varna[:-1] + ";".join([""] + ["-1"] * (longest_length - length)) + "\""
        SU.run_ct2dot(ct_dir + fname + "_unique.ct", seen_snum[sf], draw_outname_pre + ".dbn")

        # determine length of .'s needed to fill in the whole length
        OSU.system_command("sed '$s/$/&%s/' %s.dbn > %s%s_temp.dbn " % ("." * (longest_length - length), draw_outname_pre, draw_dir, varna_num))
        VIU.run_VARNA(draw_dir + str(varna_num) + "_temp.dbn", draw_outname_pre + "_structure.png", rho_varna)  # same fix this as above
        if sf == len(seen_snum) - 1 and len(seen_snum) > 1 and draw_all:  # vertical concat mult images
            v_outname = re.sub("_mult\d+", "", mult_images[0])
            VIU.vertical_image_concat(v_outname, mult_images)
            draw_outname_pre = re.findall("(.*)_structure.png$", v_outname)[0]
        if sf == len(seen_snum) - 1 or not draw_all:
            print draw_dir + str(varna_num).zfill(zero_padding) + "_structure.png"
            print "SYM LINK: " + draw_dir + str(varna_num).zfill(zero_padding) + "_structure.png"
            OSU.make_symbolic_link(re.sub("_mult\d+", "", draw_outname_pre) + "_structure.png", draw_dir + str(varna_num).zfill(zero_padding) + "_structure.png")
            VIU.convert_center_resize(draw_dir + str(varna_num).zfill(zero_padding) + "_structure.png", "1200x2800")
        OSU.remove_file(draw_dir + str(varna_num) + "_temp.dbn")
    return draw_dir + str(varna_num).zfill(zero_padding) + "_structure.png"
