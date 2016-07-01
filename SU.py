"""
Structure-related utilities (SU.py).

Version: 0.0.0
Author: Angela M Yu, 2014-2016
"""

from __future__ import division  # allows division of ints to output decimals
import os
import re
import NAU
import glob
import OSU
import math
import numpy
import shutil
import random
import multiprocessing
import collections
from itertools import repeat


def runRNAstructure_stochastic(pfsfile, ctfilename, e=1000, seed=-1, parallel=False):
    """
    Run RNAstructure executable stochastic.
    e: number of structures to sample
    seed: random seed, if not set by user, then randomly generated
    parallel: flag to use stochastic or stochastic-smp
    """
    if seed == -1:
        seed = random.randint(1, 100000000)
    if parallel:
        cmd = 'stochastic-smp %s %s -e %d --seed %d' % (pfsfile, ctfilename, e, seed)
    else:
        cmd = 'stochastic %s %s -e %d --seed %d' % (pfsfile, ctfilename, e, seed)
    cmd += " > /dev/null"
    os.system(cmd)


def runRNAstructure_partition(seqfile, pfsfile, shapefile="", constraintfile="", shape_slope=1.1, shape_intercept=-0.3, parallel=True):
    """
    Run RNAstructure executable partition.
    Requires a .seq file as input to output to pfsfile. SHAPE and constraint files are optional. SHAPE slope and intercept can be set, but are defaulted to be 1.1 and -0.3 respectively.
    parallel: flag to use partition or partition-smp
    """
    if parallel:
        cmd = 'partition-smp %s %s -si %s -sm %s' % (seqfile, pfsfile, shape_intercept, shape_slope)
    else:
        cmd = 'partition %s %s -si %s -sm %s' % (seqfile, pfsfile, shape_intercept, shape_slope)
    if shapefile != "":
        cmd += " -sh %s" % (shapefile)
    if constraintfile != "":
        cmd += " --constraint %s" % (constraintfile)
    cmd += " > /dev/null"
    os.system(cmd)


def runRNAstructure_efn2(ctfile, energyfile, parallel=False):
    """
    Run RNAstructure executable efn2.
    Requires a .ct file to output to energy file.
    parallel: flag to use efn2 or efn2-smp
    """
    if parallel:
        cmd = 'efn2-smp %s %s > /dev/null' % (ctfile, energyfile)
    else:
        cmd = 'efn2 %s %s > /dev/null' % (ctfile, energyfile)
    os.system(cmd)


def runRNAstructure_CircleCompare(predicted_ct, accepted_ct, outputfile, shape=""):
    """
    Make circle compare diagrams with RNAstructure.
    CircleCompare by default allows for flexible pairings when calculating Sensitivity and PPV, hence added the -e option for exact comparison.
    """
    cmd = 'CircleCompare -e %s %s %s' % (predicted_ct, accepted_ct, outputfile)
    if shape != "":
        cmd += " -SHAPE %s > /dev/null" % (shape)
    os.system(cmd)


def run_ct2dot(ctfile, structnum, dbnfile):
    """
    Run RNAstructure executable ct2dot which generates a .dbn file from a
    designated structure in the .ct file
    """
    structnum += 1  # ct2dot is 1-indexed and python is not
    ct2dot_cmd = 'ct2dot %s %d %s > /dev/null' % (ctfile, structnum, dbnfile)
    os.system(ct2dot_cmd)


def run_dot2ct(dbnfile, ctfile):
    """
    Run RNAstructure executable dot2ct which generates a .ct file from a
    .dbn file
    """
    dot2ct_cmd = 'dot2ct %s %s > /dev/null' % (dbnfile, ctfile)
    os.system(dot2ct_cmd)


def RNAstructure_sample(in_file_prefix, e_val, output_dir, seqfile="", shapefile="", constraintfile="", label="", num_proc=10, shape_slope=1.1, shape_intercept=-0.3, wn_tag=""):
    """
    Statistically sample one RNA, using one "mode" of RNAstructure (noshape, shape, (hard) constrained).
    Return list of structures. Uses multiple processors. Potentially can do more than one mode at a time (haven't tried it yet).
    Returns list of unique structures sampled. If a label is supplied, then it returns a list of tuples of structure and label.
    e_val: Number of total structure sampled. Will block into sets of 1000 with different random seeds.
    label: Label given to the structures obtained from this sample. Will not label if label="".
    num_proc: Number of threads to use in parallel
    """
    # Sampling mode implemented through calculation of the partition function
    if seqfile == "":
        runRNAstructure_partition(in_file_prefix+".seq", in_file_prefix+".pfs", shapefile, constraintfile, shape_slope=shape_slope, shape_intercept=shape_intercept, parallel=num_proc != 1)
    else:
        runRNAstructure_partition(seqfile, in_file_prefix+".pfs", shapefile, constraintfile, shape_slope=shape_slope, shape_intercept=shape_intercept, parallel=num_proc != 1)
    # parallelization
    e_list = [1000 if (ei+1)*1000 <= e_val else max(0, e_val - ei*1000) for ei in range(int(math.ceil(e_val/1000) + 1)) if max(0, e_val - ei*1000) > 0]
    seed_list = random.sample(xrange(1, 10000000), len(e_list))  # possible to sample from stochastic up to 100000000
    args_pool = zip(range(len(e_list)), [in_file_prefix]*len(e_list), [output_dir]*len(e_list), e_list, seed_list, repeat(wn_tag))
    structs = collections.defaultdict(int)  # holds counts for each structure
    if num_proc == 1:
        # Do no parallelization here because only 1 thread is allowed
        for i in range(0, len(args_pool)):
            tmp_structs = RNAstructure_sample_process_helper(args_pool[i])
            # Compressing into unique structures sampled
            for s in tmp_structs:
                structs[s] += 1
    else:
        for i in xrange(0, len(args_pool), num_proc):
            args_slice = args_pool[i:(i+num_proc)]
            pool = multiprocessing.Pool(processes=len(args_slice))
            p_result = pool.imap(RNAstructure_sample_process_helper, args_slice)
            for i in range(num_proc):
                tmp_structs = p_result.next(5000)
                # update structs
                for s in tmp_structs:
                    structs[s] += 1
    structs_labels = structs.keys()
    print "Unique structs sampled: " + str(len(structs_labels))
    if label != "":
        structs_labels = [(s, label) for s in structs_labels]
    return structs, structs_labels


def RNAstructure_sample_process(worker_num, in_file_prefix, output_dir, e, seed, wn_tag=""):
    """
    Process used in RNAstructure_sample. Called from RNAstructure_sample_process_helper.
    """
    wn = str(worker_num) + wn_tag
    print "Worker num: " + wn
    runRNAstructure_stochastic(in_file_prefix+".pfs", output_dir+wn+"temp.ct", e=e, seed=seed, parallel=False)
    structs = get_ct_structs(output_dir+wn+"temp.ct")
    structs_str = [",".join(s) for s in structs]
    OSU.remove_file(output_dir+wn+"temp.ct")
    return structs_str


def RNAstructure_sample_process_helper(args):
    """
    Helper to unpack the arguments and calls RNAstructure_sample_process.
    This function is needed because you can't pass multiple arguments directly to a function from Pool
    """
    return RNAstructure_sample_process(*args)


def merge_labels(list_sl, to_string=True):
    """
    Merges labels of a list of tuples where the second element in the tuple is the label.
    to_String: flag if True, then the labels are turned into a string. Else, labels are kept as a list.
    """
    sampled_structs_dict = {}
    for e in list_sl:
        struct_string = e[0]
        if isinstance(e[1], basestring):
            labels = [e[1]]
        else:
            labels = list(set(OSU.flatten_list([b.split(",") for b in e[1]])))
        if to_string:
            struct_string = ",".join(e[0])
        for l in labels:
            if struct_string not in sampled_structs_dict:
                sampled_structs_dict[struct_string] = [l]
            elif l not in sampled_structs_dict[struct_string]:
                sampled_structs_dict[struct_string].append(l)
    return list([(k, ",".join(sampled_structs_dict[k])) for k in sampled_structs_dict])


def parse_reactivity_rho(file, adapterseq, outputfile, endcut=0):
    """
    Parses reactivity file and outputs .theta and .seq, stripping off the
    adapter sequence (iteratively shrinking). Returns (positions, thetas,
    nt_sequence).
    JBL - update documentation
    """
    try:
        with open(file, 'r') as f:
            f.readline()  # throw out header
            treated_sum, untreated_sum = [[int(s)] for s in re.split("\s+", f.readline())[4:6]]  # throw out nt = 0 except for rc
            seq = []
            theta = []
            theta_cut = []
            pos = []
            adapterseq = NAU.format_rna_string(adapterseq)

            for l in f:  # parse through file
                vars = re.split('\s+', l)
                seq.append(NAU.format_rna_string(vars[3]))
                theta.append(float(vars[7]))
                pos.append(vars[2])
                untreated_sum.append(int(vars[5]))
                treated_sum.append(int(vars[4]))

            seqstring = "".join(seq)
            (seq_cut, adapter_len) = NAU.end_match_strip(seqstring, adapterseq)

            rho = calc_rho_from_theta_list(recalc_thetas(theta, 0, -adapter_len))
            # Also, endcut can be used to represent the polymerase footprint and remove these
            # bases from being used in the subsequent calculations.
            if endcut < 0:
                seq_cut = seq_cut[:endcut]
                adapter_len -= endcut
            pos = pos[:-adapter_len]
            untreated_sum = untreated_sum[:-adapter_len]
            treated_sum = treated_sum[:-adapter_len]

            rc_sum = sum(treated_sum) + sum(untreated_sum)
            rc_flag = 0 if (sum(treated_sum) + sum(untreated_sum) < 2000) else 1

            theta_cut = [str(t) for t in recalc_thetas(theta, 0, -adapter_len)]
            rho_cut = calc_rho_from_theta_list(theta_cut)
            try:
                with open(outputfile+".theta", 'w') as out:
                    out.write("\n".join(["\t".join(z) for z in zip(pos, theta_cut)]))
            except EnvironmentError:
                print "Error opening output .theta file: " + outputfile + ".theta"

            try:
                with open(outputfile+".rho", 'w') as out:
                    out.write("\n".join(["\t".join([str(zi), str(zr)]) for zi, zr in zip(pos, rho_cut)]))
            except EnvironmentError:
                print "Error opening output .rho file: " + outputfile + ".rho"

            try:
                NAU.make_seq("".join(seq_cut), outputfile+".seq")
            except EnvironmentError:
                print "Error opening output .seq file: " + outputfile + ".seq"

            return (pos, rho, theta_cut, rho_cut, seq_cut, rc_flag, rc_sum, sum(untreated_sum), sum(treated_sum))

    except EnvironmentError:
        print "Error opening reactivities file: " + file


def reactivities_to_rho_file(input_dir, adapterseq, output_dir, recalc=0, rm_temp=True, min_len=0, max_len=0):
    """
    Takes a directory of reactivities files and the adapter sequence and
    creates .theta, .rho, and .seq files. Also outputs rho_table.txt
    which is a tab delimited file of rhos. Be aware that there is numerical
    error from using floats that are present in the output. There is an option (rm_temp)
    if you do not want to keep the .theta, .rho, and .seq files.
    """
    infiles = glob.glob(input_dir + "/*_reactivities.txt")
    OSU.create_directory(output_dir)
    outname = output_dir + "rho"

    rhos = {}
    for f in infiles:
        fname = re.findall("([^/]+).txt$", f)
        output_file_prefix = output_dir+"/"+fname[0]
        pos, rho_full, theta, rho_cut, seq, rc_flag, rc_sum, ut_sum, t_sum = parse_reactivity_rho(f, adapterseq, output_file_prefix, recalc)
        if sum([float(t) for t in theta]) == 0:
            print("Not enough alignments in length {0} to calculate theta/rho. Line left blank.".format(len(pos)))
            rhos[len(pos)] = "-1\t"*len(pos) + "\n"
        else:
            rhos[len(pos)] = "\t".join([str(t) for t in rho_cut]) + "\n"

    outname += "_{0}min".format(min_len) if min_len != 0 else ""
    outname += "_{0}max".format(max_len) if max_len != 0 else ""
    outname += "_table.txt"
    with open(outname, 'w') as f:
        for key in sorted(rhos):
            if min_len == 0 or key >= min_len:  # Only apply logic if a value was supplied, but use bounds if so
                if max_len == 0 or key <= max_len:
                    f.write(rhos[key])
    if rm_temp:
        OSU.remove_file(input_dir + outname.split("/")[-1])
        shutil.move(outname, input_dir)
        OSU.remove_files_with_ext(output_dir, ".theta")
        OSU.remove_files_with_ext(output_dir, ".rho")
        OSU.remove_files_with_ext(output_dir, ".seq")
        OSU.remove_files_with_ext(output_dir, ".txt")
        os.rmdir(output_dir)
        return input_dir + os.path.basename(outname)
    return outname


def reactivities_to_reads_files(input_dir, adapterseq, min_len=0, max_len=0):
    """
    Takes a directory of reactivities files and the adapter sequence and
    creates tab delimited reads files:
    treated_mods_reads_table.txt, untreated_mods_reads_table.txt
    """
    infiles = glob.glob(input_dir + "/*_reactivities.txt")
    outname_pos = input_dir + "treated_mods_reads"
    outname_neg = input_dir + "untreated_mods_reads"
    outname_pos += "_{0}min".format(min_len) if min_len != 0 else ""
    outname_neg += "_{0}min".format(min_len) if min_len != 0 else ""
    outname_pos += "_{0}max".format(max_len) if max_len != 0 else ""
    outname_neg += "_{0}max".format(max_len) if max_len != 0 else ""
    outname_pos += "_table.txt"
    outname_neg += "_table.txt"

    reads = {}
    for file in infiles:
        with open(file, 'r') as f:
            f.readline()
            lines = f.readlines()
            seq = "".join([a.split()[3] for a in lines[1:]])
            (seq_cut, adapter_len) = NAU.end_match_strip(seq, adapterseq)
            if min_len == 0 or len(seq_cut) >= min_len:  # Only apply logic if a value was supplied, but use bounds if so
                if max_len == 0 or len(seq_cut) <= max_len:
                    reads[len(seq_cut)] = [a.split()[4:6] for a in lines[:-adapter_len]]
    with open(outname_pos, 'w') as f:
        with open(outname_neg, 'w') as f2:
            for key in sorted(reads):
                f.write("\t".join([a[0] for a in reads[key]]) + "\n")
                f2.write("\t".join([a[1] for a in reads[key]]) + "\n")

    return (outname_pos, outname_neg)


def calc_rho_from_theta(input, outputfile):
    with open(input, 'r') as f:
        theta = []
        pos = []

        for l in f:  # parse through file
            vars = re.split('\s+', l)
            theta.append(float(vars[1]))
            pos.append(vars[0])

        with open(outputfile+".rho", 'w') as out:
            rho = [repr(t * len(theta)) for t in theta]
            out.write("\n".join(["\t".join(z) for z in zip(pos, rho)]))
            return rho


def calc_rho_from_theta_list(theta):
    """
    Calculate rho reactivities from theta reactivities.
    rho = theta * length
    """
    rho = [float(t) * len(theta) for t in theta]
    return rho


def recalc_rhos(rhos, start_i, end_i):
    """
    Recalculates rho values for a subsequence of the rhos.
    It assumes that the given sequence of rhos is already normalized.
    """
    cut_rhos = [float(r) for r in rhos[start_i:end_i]]
    cut_rhos_sum = sum([float(r) for r in cut_rhos])
    cut_rhos_l = len(cut_rhos)
    return [(r*cut_rhos_l)/(cut_rhos_sum) for r in cut_rhos]


def recalc_thetas(thetas, start_i, end_i):
    """
    Recalculates theta values for a subsequence of the thetas.
    It assumes that the given sequence of rhos is already normalized.
    """
    cut_thetas = [float(r) for r in thetas[start_i:end_i]]
    cut_thetas_sum = sum([float(r) for r in cut_thetas])
    if cut_thetas_sum == 0:
        return cut_thetas
    recalced_thetas = [r/cut_thetas_sum for r in cut_thetas]
    return recalced_thetas


def rhos_list_to_file(rhos, filename):
    """ Takes a list of rhos and makes it into a filename """
    rhos = [str(r) for r in rhos]
    nts = [str(nt) for nt in range(1, len(rhos)+1)]
    lines = zip(nts, rhos)
    with open(filename, 'w') as f:
        f.write("\n".join(["\t".join(l) for l in lines]) + "\n")


def ct_list_to_file(ct, seq, filename):
    """
    Take a ct in list form and make a .ct file. The output will not
    be formatted like RNAstructure's output, but will be properly interpreted.
    """
    header = str(len(ct)) + "  " + filename
    ranges = [str(r) for r in range(len(ct)+2)]
    lines = zip(ranges[1:-1], seq, ranges[:-2], ranges[2:], ct, ranges[1:-1])
    with open(filename, 'w') as f:
        f.write(header + "\n")
        f.write("\n".join(["\t".join(l) for l in lines]) + "\n")


def cts_to_file(cts, seq, filename):
    """
    Take multiple ct's in list form and makes one .ct file. The output will not
    be formatted like RNAstructure's output, but will be properly interpreted.
    """
    #JBLQ - reminder to check this. This would be a good candidate for a unit test that
    # took an input ct, convertedit to a list, then used this to get back to a ct file.
    # could verify gives same energy with efn2 as well.
    open(filename, 'w').close()
    for i in range(len(cts)):  # loop through each ct
        ct = cts[i]
        header = "%s  %s_%s" % (str(len(ct)), filename, str(i))
        ranges = [str(r) for r in range(len(ct)+2)]
        lines = zip(ranges[1:-1], seq, ranges[:-2], ranges[2:], ct, ranges[1:-1])
        with open(filename, 'a') as f:
            f.write(header + "\n")
            f.write("\n".join(["\t".join(l) for l in lines]) + "\n")


def make_constraint_file(confile, XA, XB, XC, XD1, XD2, XE, XF1, XF2):
    """
    Writes out a contstraint file to confile, based on list input for
    all possible constraints given in list format.
    """
    XA.extend(["-1", "SS:"])
    XB.extend(["-1", "Mod:"])
    XC.extend(["-1", "Pairs:"])
    XD1_XD2 = [str(a)+" "+str(b) for a, b in zip(XD1, XD2)]
    XD1_XD2.extend(["-1 -1", "FMN:"])
    XE.extend(["-1", "Forbids:"])
    XF1_XF2 = [str(a)+" "+str(b) for a, b in zip(XF1, XF2)]
    XF1_XF2.append("-1 -1")

    with open(confile, 'w') as f:
        f.write("DS:\n")
        f.write("\n".join([str(x) for x in XA]) + "\n")
        f.write("\n".join([str(x) for x in XB]) + "\n")
        f.write("\n".join([str(x) for x in XC]) + "\n")
        f.write("\n".join(XD1_XD2) + "\n")
        f.write("\n".join([str(x) for x in XE]) + "\n")
        f.write("\n".join(XF1_XF2) + "\n")


def get_free_energy_efn2(efn2file):
    """ Returns all free energies from a ct file"""
    energy = []
    with open(efn2file, 'r') as f:
        for line in f:
            m = re.search("Energy = (.*)$", line)
            if m:
                energy.append(float(m.group(1)))
    return energy


def get_ct_structs(ctfile):
    """ Returns list of structures in a ct file in nt pairing format """
    with open(ctfile, 'r') as ct:
        structs = []
        stcurr = []
        for line in ct:
            vars = re.split('\s+', line.strip())
            if len(vars) == 6 and "ENERGY" not in vars:
                stcurr.append(vars[4])
            elif len(stcurr) > 0:
                structs.append(stcurr)
                stcurr = []
        structs.append(stcurr)  # last structure found not handled in loop
        return structs


def ct_file_to_struct_file(ctfile, outfile):
    """
    Writes a tab delimited file of structures found in the ctfile.
    """
    cts = get_ct_structs(ctfile)
    with open(outfile, 'w') as f:
        for ct in cts:
            f.write("\t".join(ct) + "\n")


def invert_scale_rho_vec(arr_rho, max_rho=-1):
    '''
    Scaling and inverting rhos based on max_rho. First projecting onto [0,1] by scaling by max_rho.
    Then inverting with 1 - so that 0s reflect un-paired nts (high rhos) and 1s represent paired nts
    (low rhos).
    '''
    if isinstance(arr_rho[0], str):
        arr_rho = [float(r) for r in arr_rho]
    if max_rho == -1:
        max_rho = max(arr_rho)
    return [0 if r >= max_rho else 1 - r/max_rho for r in arr_rho]


def scale_vec_avg1(arr, max_r=-1):
    """
    Scales a vector such that the vector's mean is 1.
    """
    arr = [float(r) for r in arr]
    if max_r != -1:
        arr = [max_r if r > max_r else r for r in arr]
    avg_v = sum(arr) / len(arr)
    if avg_v == 0:  # encountered an array of all 0's
        return [1.0 for r in arr]
    return [r/avg_v for r in arr]


def get_indices_rho_gt_c(arr_rho, max_rho, one_index=False):
    """
    Returns indices of rhos greater than a given max_rho
    """
    if one_index:
        return [i+1 for i, r in enumerate(arr_rho) if float(r) > max_rho]
    else:
        return [i for i, r in enumerate(arr_rho) if float(r) > max_rho]


def ct_struct_to_binary_vec(ct):
    """
    .ct structure(s) in numeric (string) form to 0/1 representation. 0 = unpaired, 1 = paired.
    """
    if isinstance(ct[0], int) or isinstance(ct[0], str):
        return [1 if int(a) >= 1 else 0 for a in ct]
    elif isinstance(ct[0], list):
        return [[1 if int(c) >= 1 else 0 for c in ctl] for ctl in ct]


def ct_struct_to_binary_mat(ct):
    """ .ct structure(s) in numeric or string form to 0/1 representation """
    #JBLQ - good candidate for a unit test here
    ret = []
    if isinstance(ct[0], int) or isinstance(ct[0], str):
        ret = [([0] * len(ct)) for r in xrange(len(ct))]
        ct = [int(c) for c in ct]
        for i in range(len(ct)):
            if ct[i] > 0:
                ret[i][int(ct[i])-1] = 1
    elif isinstance(ct[0], list):
        for ctl in ct:
            retl = [([0] * len(ctl)) for r in xrange(len(ctl))]
            ctl = [int(c) for c in ctl]
            for i in range(len(ctl)):
                if ctl[i] > 0:
                    retl[i][int(ctl[i])-1] = 1
            ret.append(retl)
    return ret


def cap_rho_or_ct_list(arr, max_val=-1):
    """
    Defined this to avoid using many lambda functions of this same function.
    """
    if max_val != -1:
        return [min(max_val, i) for i in arr]
    return arr


def calc_bp_distance_vector_weighted(struct, rhos, scaling_func="none", max_r=-1, invert_struct=False, paired_weight=0.5):
    """
    Calculate distance between a structure and rhos only at paired portions of the structure, normalized by number of paired and unpaired nt's in structure
    """
    if (len(struct) != len(rhos)):
        print struct
        print rhos
        raise Exception("calc_bp_distance_vector_weighted: length of lists not equal")
    struct_orig = [float(s) for s in struct]
    num_paired = sum([1 if s > 0 else 0 for s in struct_orig])
    num_unpaired = len(struct) - num_paired
    if invert_struct:  # Assumes second array is a ct file
        struct = [1 - float(a) for a in struct]
    else:
        struct = [float(a) for a in struct]
    if scaling_func == "D" or scaling_func == "K":
        scaling_func = "none"
    scaling_fns = {"U": scale_vec_avg1, "none": cap_rho_or_ct_list}
    struct = scaling_fns[scaling_func](struct, max_r)
    diffs_paired = [abs(a1 - float(a2))/num_paired for a1, a2, a3 in zip(struct, rhos, struct_orig) if a3 == 1]
    diffs_unpaired = [abs(a1 - float(a2))/num_unpaired for a1, a2, a3 in zip(struct, rhos, struct_orig) if a3 == 0]

    struct = struct_orig[:]  # copy original struct back
    return float(paired_weight)*sum(diffs_paired) + (1 - float(paired_weight))*sum(diffs_unpaired)


def calc_bp_distance_matrix(react_mat, ct_mat, endoff=0):
    """
    Calculate distance between two matrices. react_mat is assumed to have been
    made from reactivities data and elements within are scaled to fit in [0, 1]
    ct_mat is assumed to have been made from a ct file and elements are either
    0 or 1. Assumes no pseudoknots.
    endoff is used to disregard the last x bases from this distance calculation.
    For now, endoff should be negative. Ex. disregard last base => endoff=-1
    Matrices' dimensions will be changed accordingly.
    """
    if (len(react_mat) != len(ct_mat) or len(react_mat[0]) != len(ct_mat[0])):
        raise Exception("calc_bp_distance_matrix: matrix dimensions not equal")

    react_mat = numpy.matrix(react_mat)
    ct_mat = numpy.matrix(ct_mat)

    if endoff != 0:
        react_mat = react_mat[:endoff, :endoff]
        ct_mat = ct_mat[:endoff, :endoff]

    ct_mat_ut = numpy.triu(ct_mat, 1)
    react_mat_ut = numpy.triu(react_mat, 1)
    diff = numpy.absolute(numpy.subtract(ct_mat_ut, react_mat_ut))

    return numpy.sum(diff)


def calc_distances_bt_matrices(struct_matrices, endoff=0):
    """
    Returns a distance matrix between structures in struct_matrices using the matrix-based base-pair distance metric.
    """
    triu_i = numpy.triu_indices(len(struct_matrices), 1)
    distance_matrix = numpy.zeros((len(struct_matrices), len(struct_matrices)))
    for i in range(len(triu_i[0])):
        ind = (triu_i[0][i], triu_i[1][i])
        distance_matrix[ind] = calc_bp_distance_matrix(struct_matrices[ind[0]], struct_matrices[ind[1]], endoff)
        distance_matrix[ind[::-1]] = distance_matrix[ind]
    del triu_i
    return distance_matrix


def calc_benchmark_statistics_matrix(react_mat, ct_mat):
    """
    Calculate sensitivity based on matrix representation of a structure
    """
    bm_stats = {}
    react_mat = numpy.matrix(react_mat)
    ct_mat = numpy.matrix(ct_mat)
    ct_mat_ut = numpy.triu(ct_mat, 1)
    react_mat_ut = numpy.triu(react_mat, 1)
    diff = numpy.subtract(react_mat_ut, ct_mat_ut)
    FP = numpy.where(diff == 1)[0].shape[0]
    FN = numpy.where(diff == -1)[0].shape[0]
    TP = numpy.where(ct_mat_ut == 1)[0].shape[0] - FN
    print "TP: " + str(TP)
    print "FN: " + str(FN)
    print "FP: " + str(FP)
    bm_stats["F_score"] = 2*TP / float(2*TP + FN + FP)
    bm_stats["Sensitivity"] = TP / float(TP + FN)
    bm_stats["PPV"] = TP / float(TP + FP)
    print "Benchmark statistics: " + str(bm_stats)
    return bm_stats
