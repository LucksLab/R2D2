"""
Takes .pairs files and draws VARNA picture of base pairs that occur > 50% of the time

Options:
--pairs_dir       = input directory containing .pairs files
--rhos_dir        = input directory containing .rhos files
-o                = output directory

Version: 0.0.1
Author: Angela M Yu, 2014-2017

Copyright (C) 2017  Julius B. Lucks and Angela M Yu.
All rights reserved.
Distributed under the terms of the GNU General Public License, see 'LICENSE'.
"""

import SU
import OSU
import VIU
import LucksLabUtils_config
import glob
import math

# setup environment variables
LucksLabUtils_config.config("Quest_R2D2")
opts = OSU.getopts("o:", ["pairs_dir=", "rhos_dir=", "max_length="])
print opts

output_dir = OSU.create_directory(opts['-o']) + "/"
pairs_dir = opts['--pairs_dir']
rhos_dir = opts['--rhos_dir'] if '--rhos_dir' in opts else None
max_length = int(opts['--max_length']) if '--max_length' in opts else -1

rhos = {}
if rhos_dir is not None:
    for rf in glob.glob(rhos_dir+"/*_reactivities.rho"):
        # read in each rho reactivitiy spectra
        with open(rf, "r") as f:
            rho = [line.split()[1] for line in f.readlines()]
            rhos[len(rho)] = [rho, rf]  # add in rho file here
    seq_start = min(rhos.keys())
    seq_end = max(rhos.keys())
elif max_length != -1:
    seq_end = max_length
else:
    print "WARNING: assuming max length of RNA is 9999 for output purposes"
    seq_end = 9999
zero_padding = int(math.floor(math.log10(seq_end)) + 1)

for pf in glob.glob(pairs_dir+"/*.pairs"):
    # read in each pairs file
    with open(pf, "r") as f:
        seq = f.readline().strip()
        seqi = len(seq)
        if seqi in rhos:
            rho_varna = "\"" + ";".join(rhos[seqi][0]+(["-1"]*(seq_end-seqi))) + "\""
        else:
            rho_varna = "\"" + ";".join(["-1"]*(seq_end)) + "\""
        pairs = [var.split() for var in f.readlines()]
        selected_pairs = [p[:2] for p in pairs if float(p[2]) > 0.5]
        struct = ["."] * seqi
        for sp in selected_pairs:
            struct[int(sp[0])-1] = "("
            struct[int(sp[1])-1] = ")"
        with open(output_dir + str(seqi) + "_temp.dbn", "w") as f:
            f.write(">"+str(seqi)+"\n")
            f.write(seq + "\n")
            f.write("".join(struct) + "\n")
        VIU.run_VARNA(output_dir + str(seqi) + "_temp.dbn", output_dir + str(seqi).zfill(zero_padding) + "_structure.png", rho_varna)

