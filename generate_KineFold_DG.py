"""
generate_KineFold_DG.py
Parses through KineFold's .rnm output and generates a DG.dump file for
plotting.

--KineFold_dir		= Directory containing KineFold's .rnm output
--outdir		= Output directory

Optional:
--KineFold_times	= Will run KineFold a supplied number of time. Default 0
--seq_name		= Name of sequence. Default "test"
--time			= Time allowed to fold in ms. Default 160000
--pseudoknots		= Boolean to allow pseudoknots. Default False
--entanglements		= Boolean to allow entanglements. Default False
--speed			= Speed of RNA polymerase in ms per base. Default 20
--sequence		= Sequence. Default ""

Author: Angela M Yu, 2014-2017
Version: 0.0.1

Copyright (C) 2017  Julius B. Lucks, Angela M Yu.
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
opts = OSU.getopts("", ["KineFold_dir=", "outdir=", "KineFold_times=", "seq_name=", "time=", "pseudoknots=", "entanglements=", "speed=", "sequence="])
print opts
KineFold_dir = opts["--KineFold_dir"]
outdir = opts["--outdir"]
KF_times = int(opts["--KineFold_times"]) if "--KineFold_times" in opts else 0
seq_name = opts["--seq_name"] if "--seq_name" in opts else "test"
time = opts["--time"] if "--time" in opts else 160000
pseudoknots = bool(opts["--pseudoknots"] is True) if "--pseudoknots" in opts else False
entanglements = bool(opts["--entanglements"] is True) if "--entanglementss" in opts else False
speed = opts["--speed"] if "--speed" in opts else 20
sequence = opts["--sequence"] if "--sequence" in opts else ""

# create directories
dbn_outdir = OSU.create_directory(outdir + "/dbn/")
ct_outdir = OSU.create_directory(outdir + "/ct/")
efn2_outdir = OSU.create_directory(outdir + "/efn2/")

# Run KineFold
for i in range(KF_times):
    fname = "%s/%s_%s" % (KineFold_dir, seq_name, str(i))
    reqfile = SU.generate_req_dat_file(fname, sequence, time, pseudoknots, entanglements, speed)
    SU.run_KineFold(reqfile)

# parse *rnm files and create .dump file
rnm_files = glob.glob(KineFold_dir + "/*.rnm")
for rf in rnm_files:
    # parse *rnm file
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
