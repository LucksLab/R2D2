"""
generate_MFE_movie.py
Generates a MFE movie of increasing sequence length folds given a sequence.

--seq = Sequence of RNA
--outdir = Output directory

Optional:
--seq_start = 0-indexed start of sequence to be folded
--seq_end = 0-indexed start of sequence to be folded
--thetasdir = directory containing .theta files of the full pathway to be overlaid on the MFE structures. Will not be taken into account during the folding.


Author: Angela M Yu, 2014-2016
Version: 0.0.0
"""

import LucksLabUtils_config
import VIU
import OSU

LucksLabUtils_config.config("Quest_R2D2")  # set up environment

# parse command line arguments
opts = OSU.getopts("", ["seq=", "seq_start=", "seq_end=", "outdir=", "thetasdir="])
print opts
seq = opts["--seq"]
seq_start = opts["--seq_start"] if "--seq_start" in opts else -1
seq_end = opts["--seq_end"] if "--seq_end" in opts else -1
outdir = opts["--outdir"]
thetasdir = opts["--thetasdir"] if "--thetasdir" in opts else ""

# generate MFE movie
VIU.generate_MFE_CoTrans_movie(seq, seq_start, seq_end, outdir, thetasdir)
