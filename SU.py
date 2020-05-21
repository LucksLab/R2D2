"""
Structure-related utilities (SU.py).

Version: 0.0.1
Author: Angela M Yu, 2014-2016

Copyright (C) 2016  Julius B. Lucks, Angela M Yu.
All rights reserved.
Distributed under the terms of the GNU General Public License, see 'LICENSE'.
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
from sklearn.decomposition import PCA
from sklearn.manifold import MDS
from sklearn.preprocessing import scale
from sklearn.metrics import pairwise_distances
from concurrent.futures import ThreadPoolExecutor


def runRNAstructure_fold(seqfile, ctfilename, shapefile="", m=1, shape_slope=1.1, shape_intercept=-0.3, p=-99, parallel=False):
    """
    Run RNAstructure executable Fold.
    """
    if parallel:
        Fold_cmd = 'Fold-smp '
    else:
        Fold_cmd = "Fold "
    Fold_cmd += '%s %s -m %d' % (seqfile, ctfilename, m)
    if shapefile != "":
        Fold_cmd += ' -sh %s' % (shapefile)
    if p != -99:
        Fold_cmd += ' -p %d' % (p)
    Fold_cmd += ' -si %s -sm %s > /dev/null' % (shape_intercept, shape_slope)
    os.system(Fold_cmd)


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
        cmd += " -SHAPE %s" % (shape)
    cmd += " > /dev/null"
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


def RNAstructure_sample(in_file_prefix, e_val, output_dir, seqfile="", shapefile="", constraintfile="", label="", num_proc=10, shape_slope=1.1, shape_intercept=-0.3, wn_tag="", lock=None):
    """
    Statistically sample one RNA, using one "mode" of RNAstructure (noshape, shape, (hard) constrained).
    Return list of structures. Uses multiple processors. Potentially can do more than one mode at a time (haven't tried it yet).
    Returns list of unique structures sampled. If a label is supplied, then it returns a list of tuples of structure and label.
    e_val: Number of total structure sampled. Will block into sets of 1000 with different random seeds.
    label: Label given to the structures obtained from this sample. Will not label if label="".
    num_proc: Number of threads to use in parallel
    """
    # Sampling mode implemented through calculation of the partition function
    if lock is not None:
        lock.acquire()
    if seqfile == "":
        runRNAstructure_partition(in_file_prefix+".seq", in_file_prefix+".pfs", shapefile, constraintfile, shape_slope=shape_slope, shape_intercept=shape_intercept, parallel=num_proc != 1)
    else:
        runRNAstructure_partition(seqfile, in_file_prefix+".pfs", shapefile, constraintfile, shape_slope=shape_slope, shape_intercept=shape_intercept, parallel=num_proc != 1)
    if lock is not None:
        lock.release()
    # parallelization
    e_list = [1000 if (ei+1)*1000 <= e_val else max(0, e_val - ei*1000) for ei in range(int(math.ceil(e_val/1000) + 1)) if max(0, e_val - ei*1000) > 0]
    seed_list = random.sample(xrange(1, 10000000), len(e_list))  # possible to sample from stochastic up to 100000000
    args_pool = zip(range(len(e_list)), [in_file_prefix]*len(e_list), [output_dir]*len(e_list), e_list, seed_list, repeat(wn_tag), repeat(lock))
    structs = collections.defaultdict(int)  # holds counts for each structure
    
    # JBL Q - it looks like we never call this with num_proc > 1 from PCSU. Can you verify?
    # AMY - Yes, the num_proc > 1 case was from when I did not use multiprocessing.Pool in PCSU.py and wanted RNAStructure calls to be parallelized.
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


def RNAstructure_sample_process(worker_num, in_file_prefix, output_dir, e, seed, wn_tag="", lock=None):
    """
    Process used in RNAstructure_sample. Called from RNAstructure_sample_process_helper.
    """
    wn = str(worker_num) + wn_tag
    print "Worker num: " + wn
    if lock is not None:
        lock.acquire()
    runRNAstructure_stochastic(in_file_prefix+".pfs", output_dir+wn+"temp.ct", e=e, seed=seed, parallel=False)
    if lock is not None:
        lock.release()
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
    to_String: flag if True, then the structs are turned into a string. Else, structs are kept as a comma separated list.
    """
    sampled_structs_dict = {}
    for e in list_sl:
        struct_string = e[0]
        if isinstance(e[1], basestring):
            labels = [e[1]]
        else:
            labels = list(set(OSU.flatten_list([b.split(",") for b in e[1]])))
        if to_string:
            struct_string = ",".join(e[0]) #JBL Q: this is flattening the struct_string to a string, not the labels as indicated in documentation above?  #AMY: Yes, fixed it.
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
    JBL TODO - update documentation
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
            
            #cutoff adapter sequence and calculate rhos
            seqstring = "".join(seq)
            (seq_cut, adapter_len) = NAU.end_match_strip(seqstring, adapterseq)
            rho = calc_rho_from_theta_list(recalc_thetas(theta, 0, -adapter_len))

            # Remove endcut
            # For example, endcut can be used to represent the polymerase footprint and remove these
            # bases from being used in the subsequent calculations.
            if endcut < 0:
                seq_cut = seq_cut[:endcut]
                adapter_len -= endcut
            pos = pos[:-adapter_len]
            untreated_sum = untreated_sum[:-adapter_len]
            treated_sum = treated_sum[:-adapter_len]

            # rc_sum = read count sum
            # rc_flag = read count flag - used to flag lengths with low read counts
            rc_sum = sum(treated_sum) + sum(untreated_sum)
            rc_flag = 0 if (sum(treated_sum) + sum(untreated_sum) < 2000) else 1

            # Recalculate thetas and rhos with total endcut (adapter_len + endcut specified above)
            # JBL Q - why is this code block different from: 
            #      rho = calc_rho_from_theta_list(recalc_thetas(theta, 0, -adapter_len))
            # AMY - line 239 calculates rhos removing only the adapter and below removes adapter and endcut for both theta and rho calculation
            theta_cut = [str(t) for t in recalc_thetas(theta, 0, -adapter_len)]
            rho_cut = calc_rho_from_theta_list(theta_cut)
            
            # Output files
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


def binary_mat_to_binary_ct(struct_mat):
    """
    binary structure matrix (0's and 1's) to binary ct format (0's and 1's)
    """
    if not isinstance(struct_mat, numpy.matrix):
        struct_mat = numpy.matrix(struct_mat)
    return struct_mat.sum(axis=0).tolist()[0]


def binary_mat_to_ct(struct_mat):
    """
    binary structure matrix (0's and 1's) to ct format (0's for unpaired and paired nt indices (1-index))
    """
    if not isinstance(struct_mat, numpy.matrix):
        struct_mat = numpy.matrix(struct_mat)
    paired_pos = numpy.where(struct_mat == 1)
    ct = ["0"] * len(struct_mat)
    for x,y in zip(paired_pos[0], paired_pos[1]):
        ct[x] = str(y + 1)
        ct[y] = str(x + 1)
    return ct


def cap_rho_or_ct_list(arr, max_val=-1):
    """
    Defined this to avoid using many lambda functions of this same function.
    """
    if max_val != -1:
        return [min(max_val, i) for i in arr]
    return arr


def calc_bp_distance_vector_weighted(struct, rhos, scaling_func="none", max_r=-1, invert_struct=False, paired_weight=0.5):
    """
    Calculate distance between a structure and rhos, weighted by paired_weight to positions that are predicted to be
    paired and 1 - paired_weight at positions that are predicted to be unpaired.
    scaling_func = determines which distance method to use
    max_r = if not equal to -1, it caps the structure to max_r
    invert_struct = Flag to switch 0's to 1's and 1's to 0's in struct for distance calculation
    paired_weight = weight given to differences between rhos and struct at bases predicted to be paired.
    """
    if (len(struct) != len(rhos)):
        print struct
        print rhos
        raise Exception("calc_bp_distance_vector_weighted: length of lists not equal")
    struct_orig = [float(s) for s in struct]
    if len(set(struct_orig) - set([0,1])) != 0:
        raise Exception("calc_bp_distance_vector_weighted: struct must only contain 0's and 1's")
    if invert_struct:  # Assumes second array is a ct file
        struct = [1 - float(a) for a in struct]
    else:
        struct = [float(a) for a in struct]
    if scaling_func == "D" or scaling_func == "K":
        scaling_func = "none"
    scaling_fns = {"U": scale_vec_avg1, "none": cap_rho_or_ct_list}
    struct = scaling_fns[scaling_func](struct, max_r)
    diffs_paired = [abs(a1 - float(a2)) for a1, a2, a3 in zip(struct, rhos, struct_orig) if a3 == 1]
    diffs_unpaired = [abs(a1 - float(a2)) for a1, a2, a3 in zip(struct, rhos, struct_orig) if a3 == 0]

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
    if not isinstance(react_mat, numpy.matrix):
        react_mat = numpy.matrix(react_mat)
    if not isinstance(ct_mat, numpy.matrix):
        ct_mat = numpy.matrix(ct_mat)

    if react_mat.shape != ct_mat.shape:
        raise Exception("calc_bp_distance_matrix: matrix dimensions not equal")

    if endoff != 0:
        react_mat = react_mat[:endoff, :endoff]
        ct_mat = ct_mat[:endoff, :endoff]

    ct_mat_ut = numpy.triu(ct_mat, 1)
    react_mat_ut = numpy.triu(react_mat, 1)
    diff = numpy.absolute(numpy.subtract(ct_mat_ut, react_mat_ut))

    return numpy.sum(diff)


def calc_bp_distance_matrix_initializer(struct_matrices, endoff):
    """
    Used as initializer in Pool
    """
    global matrices_list, endoff_num
    matrices_list = struct_matrices
    endoff_num = endoff


def calc_bp_distance_matrix_helper_index(index):
    """
    Helper to unpack arguments and call calc_bp_distance_matrix
    args:
    0 - struct_matrix 1
    1 - struct_matrix 2
    2 - endoff
    3 - index
    """
    return (index, calc_bp_distance_matrix(matrices_list[index[0]], matrices_list[index[1]], endoff_num))


def calc_distances_bt_matrices(struct_matrices, endoff=0, n_jobs=1):
    """
    Returns a distance matrix between structures in struct_matrices using the matrix-based base-pair distance metric.
    """
    triu_i = numpy.triu_indices(len(struct_matrices), 1)
    distance_matrix = numpy.zeros((len(struct_matrices), len(struct_matrices)))
    if n_jobs == 1:
        for i in range(len(triu_i[0])):
            ind = (triu_i[0][i], triu_i[1][i])
            distance_matrix[ind] = calc_bp_distance_matrix(struct_matrices[ind[0]], struct_matrices[ind[1]], endoff)
            distance_matrix[ind[::-1]] = distance_matrix[ind]
    else:
        ind = zip(triu_i[0], triu_i[1])
        """
        pool = multiprocessing.Pool(processes=n_jobs, initializer = calc_bp_distance_matrix_initializer, initargs=(struct_matrices, endoff))
        pool_results = pool.imap_unordered(calc_bp_distance_matrix_helper_index, [i for i in ind], chunksize=10000)
        pool.close()
        pool.join()
        """
        """
        calc_bp_distance_matrix_initializer(struct_matrices, endoff)
        with ThreadPoolExecutor(max_workers=n_jobs) as executor:
            futures = executor.map(calc_bp_distance_matrix_helper_index, [i for i in ind]) # start
            pool_results = [future for future in futures] # wait for results
        """
        calc_bp_distance_matrix_initializer(struct_matrices, endoff)
        """
        # 4120765 109 nt
        for big_chunk in xrange(0, len(ind), 180000000):
            cur_ind = ind[big_chunk:(big_chunk+180000000)]
            pool = multiprocessing.Pool(processes=n_jobs)
            pool_results = pool.imap_unordered(calc_bp_distance_matrix_helper_index, [i for i in ind], chunksize=min(10000, len(ind) / (n_jobs * 4)))
            pool.close()
            pool.join()
            del pool
            for curr_ind, res in pool_results:
                distance_matrix[curr_ind] = res
                distance_matrix[curr_ind[::-1]] = res
            del pool_results, curr_ind, res
            n_jobs = max(2, n_jobs/2)  # protect memory by forking less each round, TODO: dynamically determine this
        """
        pool = multiprocessing.Pool(processes=n_jobs, maxtasksperchild=180000000)
        pool_results = pool.imap_unordered(calc_bp_distance_matrix_helper_index, [i for i in ind], chunksize=min(10000, len(ind) / (n_jobs * 4)))
        pool.close()
        pool.join()
        del pool
        for curr_ind, res in pool_results:
            distance_matrix[curr_ind] = res
            distance_matrix[curr_ind[::-1]] = res
        del pool_results, curr_ind, res
    del triu_i
    return distance_matrix


def run_PCA(X, reactivities_prefix, center=False, scale_std=False):
    """
    PCA

    WARNING: Python's implementation of PCA may not give accurate coordinates for same structures
    """
    pca = PCA(n_components=2)
    if not isinstance(X, numpy.matrix):
        X = numpy.matrix(X)
    if center is True or scale_std is True:
        output_prefix = []
        X = scale(X, with_mean=center, with_std=scale_std)
        if center is True:
            output_prefix.append("centered")
        if scale_std is True:
            output_prefix.append("scaled")
        output_prefix = "_".join(output_prefix) + "_"
    else:
        output_prefix = ""
    principal_components = pca.fit_transform(X)
    with open("%s_%sPCA_coords.txt" % (reactivities_prefix, output_prefix), "w") as f:
        f.write("\n".join(["\t".join([str(coord) for coord in row]) for row in principal_components.tolist()]) + "\n")
    return principal_components


def run_MDS_mat(X_mats, reactivities_prefix, p=1, MDS_p=1):
    """
    WARNING: MDS implementation in python may have same eigenvalue bug as cmdscale() in R
    https://stat.ethz.ch/pipermail/r-sig-ecology/2010-July/001390.html
    """
    distances = calc_distances_bt_matrices(X_mats, n_jobs=p)
    print "Finished run_MDS_mat distances, %s processes" % (p)
    model = MDS(n_components=2, dissimilarity='precomputed', n_jobs=MDS_p, random_state=1)  # don't need random_state?
    print "Starting fit_transform, %s processes" % (MDS_p)
    mds_coords = model.fit_transform(distances)
    print "Finished fit_transform"
    del distances, model
    with open("%s_MDS_mat_coords.txt" % (reactivities_prefix), "w") as f:
        f.write("\n".join(["\t".join([str(coord) for coord in row]) for row in mds_coords.tolist()]) + "\n")
    return mds_coords


def run_MDS_ct(X, reactivities_prefix, p=1, MDS_p=1):
    """
    Runs a ct version of MDS
    X - assumed to be vectors of binary ct structures in each row

    WARNING: MDS implementation in python may have same eigenvalue bug as cmdscale() in R
    https://stat.ethz.ch/pipermail/r-sig-ecology/2010-July/001390.html
    """
    if not isinstance(X, numpy.matrix):
        X = numpy.matrix(X)
    print "Starting run_MDS_ct distances, %s processes" % (p)
    distances = pairwise_distances(X, metric='manhattan', n_jobs=p)
    print "Finished run_MDS_ct distances"
    model = MDS(n_components=2, dissimilarity='precomputed', n_jobs=MDS_p, random_state=1)
    print "Starting fit_transform, %s processes" % (MDS_p)
    mds_coords = model.fit_transform(distances)
    print "Finished fit_transform"
    del distances, model
    with open("%s_MDS_ct_coords.txt" % (reactivities_prefix), "w") as f:
        f.write("\n".join(["\t".join([str(coord) for coord in row]) for row in mds_coords.tolist()]) + "\n")
    return mds_coords


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
    TN_mat = sum(range(ct_mat_ut.shape[0])) - FP - FN - TP
    TN_ct = ct_mat_ut.shape[0] - FP - FN - TP
    print "TP: " + str(TP)
    print "FN: " + str(FN)
    print "FP: " + str(FP)
    print "TN_mat: " + str(TN_mat)
    print "TN_ct: " + str(TN_ct)
    bm_stats["F_score"] = 2*TP / float(2*TP + FN + FP) if TP + FN + FP != 0 else float('nan')
    bm_stats["Sensitivity"] = TP / float(TP + FN) if TP + FN != 0 else float('nan')
    bm_stats["PPV"] = TP / float(TP + FP) if TP + FP != 0 else float('nan')
    bm_stats["TP"] = TP
    bm_stats["FN"] = FN
    bm_stats["FP"] = FP
    bm_stats["TN_mat"] = TN_mat
    bm_stats["TN_ct"] = TN_ct
    print "Benchmark statistics: " + str(bm_stats)
    return bm_stats


def run_KineFold(reqfile):
    """
    Calls KineFold with supplied reqfile.
    """
    OSU.system_command("kinefold_long_static %s -noprint" % (reqfile))


def generate_req_dat_file(fpre, sequence, time=160000, pseudoknots=False, entanglements=False, speed=20, seed=-1):
    """
    Generates a .req and .dat file to run KineFold. Must provide file name prefix (fpre) and sequence.
    Optionally can supply time allowed to fold, allow pseudoknots, allow entanglements,RNA polymerase
    speed, and seed.
    """
    extensions = [".p", ".e", ".rnm", ".rnms", ".rnml", ".rnm2", ".dat"]
    # generate random seed if neeeded
    if seed == -1:
        seed = random.randint(1, 100000000)

    # create .req file
    with open(fpre + ".req", "w") as f:
        f.write("%s\t\t# randomseed\n" % (seed))
        f.write("\n".join(["%s%s" % (fpre, ext) for ext in extensions]) + "\n")
        f.write("0\t\t# 0=RNA ; 1=DNA\n")
        f.write("6.3460741\t# helix minimum free energy in kcal/mol: 6.3460741=10kT\n")
        f.write("10000000\t# NA\n")
        f.write("%s\t\t# folding time requested in msec\n" % (time))
        f.write("%s\t\t# pseudoknots   1=yes 0=no\n" % (int(pseudoknots)))
        f.write("%s\t\t# entanglements 1=yes 0=no\n" % (int(entanglements)))
        f.write("2 %s\t\t# simulation type: 1=renaturation; 2 20=cotrans. @ 20msec/nt\n" % (speed))
        f.write("\t\t# add T i j k or F i j k options here\n\n")
        f.write("%s\n%s.zip\n<SEQNAME>%s\n<BASE>%s\n" % (fpre, fpre, fpre, fpre))
        f.write("<SEQUENCE>%s\n<ZIPFILE>%s.zip\n" % (sequence, fpre))

    # makes .dat file
    with open(fpre + ".dat", "w") as f:
        f.write("< %s\n%s\n" % (fpre, sequence))

    return fpre + ".req"


def get_rnm_structs_dbn(rnmfile, outputdir, return_time=False, total_time=-1):
    """
    Takes the .rnm output from KineFold and creates associated .dbn files in
    the output directory. Outputs the names of the output files and respective
    KineFold free energies in order.
    If return_time is set, then it will return the time spent in each
    structure in ms. total_time must also be set in ms, ex. 40000 for 40 sec. 
    """
    with open(rnmfile, 'r') as f:
        fpre = re.findall("([^/]+).rnm$", rnmfile)[0]
        seq = ""
        seqr = ""
        length = "0"
        count = 1
        files = []
        energy_path = []
        if return_time:
            # KineFold seems to start at a longer length than 0 and from ssRNA
            prev_time = -1
            time_spent = []
        for line in f:
            seq_tmp = re.match("^([ACUG\s\[\]\^]+)\| ([-\d\.]+).*after ([\d\.]+).* (\d+) over", line)
            h_rep = re.match("^([ \-\d\']+) H\s+Helix numbering$", line)
            if seq_tmp:
                seqr = seq_tmp.group(1)
                energy_path.append(seq_tmp.group(2))
                seq = re.sub("(\\]|\\[|\^)", ' ', seqr)
                seq = "".join(seq.split())
                if length == str(len(seq)):
                    count += 1
                else:
                    count = 1
                    length = str(len(seq))
                if return_time:
                    curr_time = float(seq_tmp.group(3))
                    if prev_time >= 0:
                        time_spent.append(curr_time - prev_time)
                    prev_time = curr_time
            if h_rep:
                db = rnm_to_dotbracket(seqr, h_rep.group(1))
                fname = "_".join([fpre, length, str(count) + ".dbn"])
                with open(outputdir + "/" + fname, "w") as w:
                    w.write("> " + fname + "\n")
                    w.write(seq + "\n")
                    w.write(db + "\n")
                files.append(outputdir + "/" + fname)
        if return_time:
            # add in last time and check for potential user input error
            assert total_time >= prev_time, "Total simulation time is less than time found in .rnm file"
            time_spent.append(total_time - prev_time)
            return files, energy_path, time_spent
        else:
            return files, energy_path


def rnm_to_dotbracket(seql, h_rep):
    """
    Takes the parsed two lines in an .rnm file and creates the associated dot
    bracket notation.
    """
    seql = seql.strip()
    h_rep = h_rep.strip()
    seql = re.sub('\^', '][', seql)
    seql = re.sub("(\\]|\\[)", r' \1 ', seql)
    helices_seq = re.findall("\\[[\sACUG]+\\]", seql)
    digit_re = re.compile("\d")
    left_par_re = re.compile("\'")
    h_rep = [h for h in h_rep.split(" ") if bool(digit_re.search(h))]
    for s in range(len(helices_seq)):
        if bool(left_par_re.search(h_rep[s])):
            seql = seql.replace(helices_seq[s], ')' * ((len(helices_seq[s])//2) - 1), 1)
        else:
            seql = seql.replace(helices_seq[s], '(' * ((len(helices_seq[s])//2) - 1), 1)
    seql = "".join(seql.split())
    seql = re.sub("[ACUG]", ".", seql)
    return seql
