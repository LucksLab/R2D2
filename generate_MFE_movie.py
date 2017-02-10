"""
generate_MFE_movie.py
Generates a MFE movie of increasing sequence length folds given a sequence.

--seq = Sequence of RNA
--outdir = Output directory

Optional:
--seq_start 	= 0-indexed start of sequence to be folded
--seq_end 	= 0-indexed start of sequence to be folded
--rhos_dir 	= directory containing .rho files of the full pathway to be overlaid on the MFE structures. Will not be taken into account during the folding unless --SHAPE_direct == True
--SHAPE_direct	= Turn on SHAPE directed folding. Default False.


Author: Angela M Yu, 2014-2017
Version: 0.0.1

Copyright (C) 2016, 2017  Julius B. Lucks, Angela M Yu.
All rights reserved.
Distributed under the terms of the GNU General Public License, see 'LICENSE'.
"""

import LucksLabUtils_config
import VIU
import OSU

LucksLabUtils_config.config("Quest_R2D2")  # set up environment

# parse command line arguments
opts = OSU.getopts("", ["seq=", "seq_start=", "seq_end=", "outdir=", "rhos_dir=", "SHAPE_direct="])
print opts
seq = opts["--seq"]
outdir = opts["--outdir"]
seq_start = opts["--seq_start"] if "--seq_start" in opts else -1
seq_end = opts["--seq_end"] if "--seq_end" in opts else -1
rhos_dir = opts["--rhos_dir"] if "--rhos_dir" in opts else ""
SHAPE_direct = bool(opts["--SHAPE_direct"] == True) if "--SHAPE_direct" in opts else False
print outdir
# generate MFE movie
VIU.generate_MFE_CoTrans_movie(seq, outdir, seq_start, seq_end, rhos_dir, SHAPE_direct)
