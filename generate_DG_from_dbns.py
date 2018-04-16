"""
Takes .dbn files and generates *DG.dump files containing DG of the structures in the dbn files

Options:
--dbn_dir       = input directory containing .dbn files
-o              = output directory
--out_prefix    = output file prefix

Version: 0.0.1
Author: Angela M Yu, 2014-2017

Copyright (C) 2017  Julius B. Lucks and Angela M Yu.
All rights reserved.
Distributed under the terms of the GNU General Public License, see 'LICENSE'.
"""

import SU
import OSU
import LucksLabUtils_config
import glob

# setup environment variables
LucksLabUtils_config.config("Quest_R2D2")
opts = OSU.getopts("o:", ["dbn_dir=", "out_prefix="])
print opts

output_dir = OSU.create_directory(opts['-o'])
dbn_dir = opts['--dbn_dir']
out_prefix = opts['--out_prefix']

dbns = {}
for dbnf in glob.glob(dbn_dir+"/*.dbn"):
    # read in each rho reactivitiy spectra
    with open(dbnf, "r") as f:
        dbn = f.readlines()[-1].strip()  # last line with dotbracket
        SU.run_dot2ct(dbnf, output_dir + "temp.ct")
        SU.runRNAstructure_efn2(output_dir + "temp.ct", output_dir + "temp.efn2", parallel=False)
        dg = SU.get_free_energy_efn2(output_dir + "temp.efn2")[0]
        dbns[len(dbn)] = dg  # save in dict, nt : DG

OSU.remove_files([output_dir + "temp.ct", output_dir + "temp.efn2"])

with open(output_dir + out_prefix + "_DG.dump", "w") as f:
    f.write("\t".join(["nt", "DG", "consensus"]) + "\n")
    for key in sorted(dbns):
        f.write("\t".join([str(key), str(dbns[key]), "1"]) + "\n")
