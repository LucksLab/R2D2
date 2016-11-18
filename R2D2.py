"""
R2D2 (Reconstructing RNA Dynamics with Data)

Uses cotranscriptional SHAPE-Seq reactivities to determine the cotranscriptional folding pathway of an RNA.

Version: 0.0.1
Author: Angela M Yu, 2014-2016

Copyright (C) 2016  Julius B. Lucks, Angela M Yu.
All rights reserved.
Distributed under the terms of the GNU General Public License, see 'LICENSE'.
"""

import PCSU
import OSU
import VIU
from multiprocessing import Pool
from collections import defaultdict
import glob
from itertools import repeat
import math


class R2D2:
    """
    R2D2 (Reconstructing RNA Dynamics with Data)

    Uses Cotranscriptional SHAPE-Seq reactivities to determine the
    cotranscriptional folding pathway of an RNA that is supported by the data.
    """

    def __init__(self, inputdir, outputdir, adapterseq,
                 shape_slope=1.1, shape_intercept=-0.3, p=1, e=1000, endcut=0,
                 constrained_c=2.6, scale_rho_max=2.3, draw_all=True, most_count_tie_break=True,
                 weight_paired=0.3, scaling_func="K", cap_rhos=True, pol_fp=0):
        """
        R2D2 (Reconstructing RNA Dynamics with Data) class definition.

        OPTIONS:
        # inputdir 	= Input directory containing reactitivy files
        # outputdir 	= Output directory
        # adapterseq 	= Adapter sequence
        # shape_slope 	= Slope used with SHAPE restraints in RNAstructure
        # shape_intercept = Intercept used with SHAPE restraints in RNAstructure
        # p 		= Number of threads to use
        # e 		= Number of structures sampled by each sampling method
        # endcut 	= Number of bases to be ommitted from calculations from the 3' end of the length
                            Ex. Ignoring the last base would correspond to endcut=-1
        # pol_fp 	= Number of bases occupied by polymerase. Globally sets for each length.
                            Ex. Polymerase foot print of size 9 => pol_fp=9
        # constrained_c = Any rho value greater than or equal to this value will be forced as unpaired when sampling with hard constraints
        # scale_rho_max = If True or the 'D' distance is used, this value is the max value and all values greater than it are set to this max value
        # draw_all 	= Whether or not to show all best structures in the movie
        # most_count_tie_break = Flag to use count of structures as the tie breaking criteria when showing only one structure in the movie (draw_all=False). If True, then the structure sampled the most from the pool of sampled structures is drawn.
        # weight_paired = Weight given to paired regions in distance calculations.
        # scaling_func 	= Choice of distance function when choosing the best structure:
                            D: Bound to be between [0,1]
                            U: Rescale sampled structures to average to 1
                            K: Keep sampled structures and reactivities values. If cap_rhos is True, then reactivities will be capped.
        # cap_rhos 	= Flag to have a max cutoff for reactivities when calculating distances for choosing the best structure
        """
        self.file_data = defaultdict(dict)
        self.shape_slope = shape_slope
        self.shape_intercept = shape_intercept
        self.output_dir = outputdir
        self.input_dir = inputdir
        self.adapterseq = adapterseq
        self.p = p
        self.e = e
        self.endcut = endcut
        self.pol_fp = pol_fp
        self.constrained_c = constrained_c
        self.scale_rho_max = scale_rho_max
        self.draw_all = draw_all
        self.most_count_tie_break = most_count_tie_break
        self.weight_paired = weight_paired
        self.struct_distances = {}
        self.scaled_rhos = {}
        self.scaling_func = scaling_func
        self.cap_rhos = cap_rhos
        self.run()

    def run(self):
        """
        The main routine of R2D2.
        Parses reactivities output from Spats and controls the calls to PCSU.
        """

        max_best_states = -1  # max number of best states across the lengths
        OSU.create_directory(self.output_dir)
        ct_dir = OSU.create_directory(self.output_dir + "/ct/")
        pickle_dir = OSU.create_directory(self.output_dir + "/pickles/")
        infiles = glob.glob(self.input_dir + "/*_reactivities.txt")

        # Pre-processing all input reactivities files - trimming adapter, recalculating thetas, calculating rhos
        max_best_states = 0
        rhos = {}
        rhos_cut = {}

        # Set up and run parallized calculations on each length
        args_pool = zip(infiles, repeat(self.output_dir), repeat(ct_dir),
                        repeat(pickle_dir), repeat(self.adapterseq),
                        repeat(self.endcut), repeat(self.pol_fp),
                        repeat(self.e), repeat(self.constrained_c),
                        repeat(self.cap_rhos), repeat(self.scale_rho_max),
                        repeat(self.scaling_func), repeat(self.weight_paired))
        print "run args_pool length: " + str(len(args_pool))

        if self.p > 1:  # start pool if multithread
            pool = Pool(processes=self.p)
            for length_key, file_data_length_key, struct_distances_length, num_min_states, rho, rho_cut in pool.imap(PCSU.run_cotrans_length_helper, args_pool):
                print "done length_key: " + str(length_key)
                if max_best_states < num_min_states:
                    max_best_states = num_min_states
                self.file_data[length_key] = file_data_length_key
                self.struct_distances[length_key] = struct_distances_length
                rhos[length_key + self.endcut] = "\t".join([str(r) for r in rho]) + "\n"
                rhos_cut[length_key + self.endcut] = "\t".join([str(r) for r in rho_cut]) + "\n"
        else:  # no multiprocessing
            for args_slice in args_pool:
                length_key, file_data_length_key, struct_distances_length, num_min_states, rho, rho_cut = PCSU.run_cotrans_length_helper(args_slice)
                print "done length_key: " + str(length_key)
                if max_best_states < num_min_states:
                    max_best_states = num_min_states
                self.file_data[length_key] = file_data_length_key
                self.struct_distances[length_key] = struct_distances_length
                rhos[length_key + self.endcut] = "\t".join([str(r) for r in rho]) + "\n"
                rhos_cut[length_key + self.endcut] = "\t".join([str(r) for r in rho_cut]) + "\n"

        # Output the rho reactivity matrix
        with open(self.output_dir + "rho_table.txt", 'w') as f:
            print "sorted(rhos): " + str(len(rhos.keys()))
            for key in sorted(rhos):
                f.write(rhos[key])
        with open(self.output_dir + "rho_table_cut.txt", 'w') as f:
            print "sorted(rhos): " + str(len(rhos_cut.keys()))
            for key in sorted(rhos_cut):
                f.write(rhos_cut[key])

        # organizing files into their respective directories
        for file_ext in ["rho", "theta", "seq", "pfs", "con", "efn2"]:
            OSU.create_directory(self.output_dir+file_ext+"_dir/")
            OSU.system_command("mv %s/*%s %s/%s_dir/" % (self.output_dir, file_ext, self.output_dir, file_ext))

        self.generate_output()  # generate majority of output

    def generate_output(self):
        """
        Majority of Cotranscriptional SHAPE-Seq output is created here. This includes the DG plot, best structure images, and the movie of the best structures.
        """
        draw_dir = OSU.create_directory(self.output_dir + "/draw/")
        nn_dir = OSU.create_directory(self.output_dir + "/nn/")
        OSU.create_directory(nn_dir + "distances/")
        sorted_lengths = sorted(self.file_data)
        zero_padding = int(math.floor(math.log10(sorted_lengths[-1])) + 1)

        # Parallelized function calls to generate DG plot, distance matrices for clustering of structures, and creating images of minimum distance structures
        draw_struct_nums = [length for length in sorted_lengths if self.file_data[length]["rc_flag"]]
        draw_args_pool = zip([self.file_data[dsn] for dsn in draw_struct_nums],
                             draw_struct_nums, repeat(sorted_lengths[-1]),
                             range(1, len(draw_struct_nums)+1), repeat(zero_padding),
                             repeat(draw_dir), repeat(self.output_dir + "/ct/"),
                             repeat(self.draw_all), repeat(self.most_count_tie_break))
        args_pool = [(PCSU.generate_DG_output, (self, 1, sorted_lengths[-1]))] + zip(repeat(PCSU.generate_best_struct_images), draw_args_pool)

        if self.p == 1:
            for i in range(len(args_pool)):
                PCSU.calculate_function_helper(args_pool[i])
        else:
            PCSU.run_output_multiprocessing_pool(PCSU.calculate_function_helper, args_pool, self.p)

        if not OSU.check_file_exists(self.output_dir+"/DG_state_plot.pdf"):  # Weird error on quest that it will ignore this command if sample size is very large
            PCSU.generate_DG_output(self, 1, sorted_lengths[-1])


        # Use ffmpeg on varna_num.png's to create video of the minimum distance folding pathway
        OSU.make_symbolic_link(draw_dir + str(len(draw_struct_nums)).zfill(zero_padding) + "_structure.png", draw_dir + str(len(draw_struct_nums)+1).zfill(zero_padding) + "_structure.png")  # ffmpeg needs a duplicate of the last frame
        VIU.generate_movie(draw_dir + "%%%dd_structure.png" % (zero_padding), self.output_dir + "/movie.mp4")
        return
