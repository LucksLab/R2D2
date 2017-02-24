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
--make_DG_dump	= Turn on generation of DG dump file from folding pathway. Default False


Author: Angela M Yu, 2014-2017
Version: 0.0.1

Copyright (C) 2016, 2017  Julius B. Lucks, Angela M Yu.
All rights reserved.
Distributed under the terms of the GNU General Public License, see 'LICENSE'.
"""

import LucksLabUtils_config
import VIU
import OSU
import SU
import re

LucksLabUtils_config.config("Quest_R2D2")  # set up environment

# parse command line arguments
opts = OSU.getopts("", ["seq=", "seq_start=", "seq_end=", "outdir=", "rhos_dir=", "SHAPE_direct=", "make_DG_dump="])
print opts
seq = opts["--seq"]
outdir = opts["--outdir"]
seq_start = int(opts["--seq_start"]) if "--seq_start" in opts else -1
seq_end = int(opts["--seq_end"]) if "--seq_end" in opts else -1
rhos_dir = opts["--rhos_dir"] if "--rhos_dir" in opts else ""
SHAPE_direct = bool(opts["--SHAPE_direct"] == "True") if "--SHAPE_direct" in opts else False
make_DG_dump = bool(opts["--make_DG_dump"] == "True") if "--make_DG_dump" in opts else False

# generate MFE movie
VIU.generate_MFE_CoTrans_movie(seq, outdir, seq_start, seq_end, rhos_dir, SHAPE_direct)

# generate DG dump file
if make_DG_dump:
    ct_dir = outdir + "/ct/"
    efn2_dir = OSU.create_directory(outdir + "/efn2/")
    name_nums = range(1, len(seq)+1)
    if seq_start != -1:
        name_nums = name_nums[seq_start - 1:]
    if seq_end != -1:
        name_nums = name_nums[:seq_end - seq_start + 2]

    # write dumpfile header
    fname_dump = outdir + re.match(".*\/(.*)$", outdir.rstrip("/")).group(1) + "_DG_state_plot.dump"
    with open(fname_dump, "w") as f:
        f.write("nt\tDG\n")

    for n in name_nums:
        ct_file = "%s%s.ct" % (ct_dir, n)
        if not OSU.check_file_exists(ct_file):
            continue
        energy_file = "%s%s.efn2" % (efn2_dir, n)
        SU.runRNAstructure_efn2("%s%s.ct" % (ct_dir, n), energy_file)
        energy = SU.get_free_energy_efn2(energy_file)[0]
        with open(fname_dump, "a") as f:
            f.write("%s\t%s\n" % (n, energy))
