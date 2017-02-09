"""
generate_KineFold_DG.py
Parses through KineFold's .rnm output and generates a DG.dump file for
plotting.

--KineFold_dir	= Directory containing KineFold's .rnm output
--outdir	= Output directory


Author: Angela M Yu, 2014-2016
Version: 0.0.1

Copyright (C) 2016  Julius B. Lucks, Angela M Yu.
All rights reserved.
Distributed under the terms of the GNU General Public License, see 'LICENSE'.
"""

import LucksLabUtils_config
import OSU
import SU
import glob
import re

LucksLabUtils_config.config("Quest_R2D2")  # set up environment

# parse command line arguments
opts = OSU.getopts("", ["KineFold_dir=", "outdir="])
print opts
KineFold_dir = opts["--KineFold_dir"]
outdir = opts["--outdir"]
dbn_outdir = OSU.create_directory(outdir + "/dbn/")
ct_outdir = OSU.create_directory(outdir + "/ct/")
efn2_outdir = OSU.create_directory(outdir + "/efn2/")

# parse *rnm files and create .dump file
rnm_files = glob.glob(KineFold_dir + "*.rnm")
for rf in rnm_files:
    # parse *rnm file
    print rf
    kf_dbns, kf_energy_path = SU.get_rnm_structs_dbn(rf, dbn_outdir)

    # create dump with header
    fname_dump = re.match(".*\/(.*).rnm$", rf).group(1) + "_DG_state_plot.dump"
    with open("%s/%s" % (outdir, fname_dump), "w") as f:
        f.write("nt\tDG\tKineFold_DG\n")

    # fill in dump file
    assert(len(kf_dbns) == len(kf_energy_path))
    for dbn, kf_energy in zip(kf_dbns, kf_energy_path):
        fpre = re.match(".*\/(.*_(\d+)_\d+).dbn$", dbn)
        fname_ct = ct_outdir + fpre.group(1) + ".ct"
        fname_efn2 = efn2_outdir + fpre.group(1) + ".efn2"

        SU.run_dot2ct(dbn, fname_ct)
        SU.runRNAstructure_efn2(fname_ct, fname_efn2)
        energy = SU.get_free_energy_efn2(fname_efn2)[0]
        with open("%s/%s" % (outdir, fname_dump), "a+") as d:
            d.write("\t".join([fpre.group(2), str(energy), kf_energy]) + "\n")

# plotting
