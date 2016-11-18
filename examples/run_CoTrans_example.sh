#!/usr/bin/bash

# Example of how to run R2D2 (Reconstructing RNA Dynamics from Data) on a cotranscriptional SHAPE-Seq data set

# Options:
#--in_dir : Input directory containing reactivities files
#--out_dir : Output directory
#--adapter : Adapter sequence
#--e : Size of sample to be used for each of the sampling methods
#--p : Number of threads allowed to use, default 1
#--endcut : Removes 3' indices based on value passed. Ex. --endcut = -1 => removes last base from input reads and reactivities
#--pol_fp : Remove 3' indices based on length of RNA polymerase footprint
#--constrained_c : Parameter for hard-constrained sampling. Any rho value greater than this value is forced to be unpaired.
#--scale_rho_max : Parameter for rescaling rhos such that rhos are capped to this value.
#--draw_all : Flag for whether or not to draw all possible best states for the best structure path video.
#--most_count_tie_break : When making the video of the best structure path, this flag determines if the structure sampled the most number of times is used instead of all best structures. This flag is only relevant if --draw_all is False.  
#--weight_paired : Weight parameter for weighted distance calculation. 
#--scaling_func = Choice of distance function when choosing the best structure:
#                  D: Bound to be between [0,1]
#                  U: Rescale sampled structures to average to 1
#                  K: Keep sampled structures and reactivities values. If cap_rhos is True, then reactivities will be capped.
#--cap_rhos = Flag to have a max cutoff when calculating distances for choosing the best structure


# Output:
#draw/ : directory with output related to making structure images and videos
#ct/ : directory of structures sampled
#pickles/ : directory of pickled data
#movie.mp4 : Video of structures along the best structure path
#pfs/ : directory of partition functions
#seq/ : directory of sequence files
#theta/ : directory of theta files
#rho/ : directory of rho files
#*dump : output to be used for plotting in R
#rho_table.txt : table of rho values sorted by length
#rho_table_cut.txt : table of rho values sorted by length after removing 3' end nucleotides specified by --endcut and --pol_fp
#./CoTrans_example_output/DG_state_plot.pdf : Plot of delta G vs length. Cotranscriptional folding pathway is denoted with red.

cwd=`pwd`

python ../analyze_cotrans_SHAPE-Seq.py --in_dir "$cwd/reactivities_3/" --out_dir "$cwd/CoTrans_example_output/" --adapter "CTGACTCGGGCACCAAGG" --e 100 --endcut 0 --constrained_c "1.8" --scale_rho_max "1.8" --draw_all "False" --most_count_tie_break "False" --weight_paired "0.5" --scaling_func "K" --cap_rhos "True" --pol_fp 14 --p 1

python ../analyze_cotrans_SHAPE-Seq.py --in_dir "$cwd/reactivities_1/" --out_dir "$cwd/CoTrans_example_output_1/" --adapter "CTGACTCGGGCACCAAGG" --e 100 --endcut 0 --constrained_c "1.8" --scale_rho_max "1.8" --draw_all "False" --most_count_tie_break "False" --weight_paired "0.5" --scaling_func "K" --cap_rhos "True" --pol_fp 14 --p 1
