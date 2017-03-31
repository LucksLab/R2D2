"""
Combines DG_state_plot.dump files. Only lengths where rc_flag == 1 are kept.
This is a script for the initial R2D2 paper in which the 100 times R2D2 were
run went in the order of "3_times_dirs", "50_times_dir", followed by
"47_times_dir". Will also unpack DG_state_plot.dump from results_except_draw.tgz
if DG_state_plot.dump not found in folder.

Options:
-o		= Output directory
--3_times_dirs	= List of first 3 directories separated by ","
--50_times_dir	= Directory which contains output directories from 50_times
--47_times_dir	= Directory which contains output directories from 47_times
--file_prefix   = Prefix of runs in 50_time_dir and 47_times_dir

Version: 0.0.1
Author: Angela M Yu, 2014-2016

Copyright (C) 2016  Julius B. Lucks and Angela M Yu.
All rights reserved.
Distributed under the terms of the GNU General Public License, see 'LICENSE'.
"""

import OSU
import LucksLabUtils_config
from collections import defaultdict

# setup environment variables
LucksLabUtils_config.config("Quest_R2D2")
opts = OSU.getopts("o:", ["3_times_dirs=", "50_times_dir=", "47_times_dir=", "100_times_dir=", "file_prefix="])
print opts

output_dir = OSU.create_directory(opts['-o'])
file_prefix = opts["--file_prefix"]

if "--100_times_dir" in opts:
    times_dirs = ["%s/%s%s/" % (opts["--100_times_dir"], file_prefix, i) for i in range(1, 101)]
elif "--3_times_dirs" in opts and "--50_times_dir" in opts and "--47_times_dir" in opts:
    times_dirs = opts["--3_times_dirs"].split(",")
    times_dirs += ["%s/%s%s/" % (opts["--50_times_dir"], file_prefix, i) for i in range(1, 51)]
    times_dirs += ["%s/%s%s/" % (opts["--47_times_dir"], file_prefix, i) for i in range(1, 48)]
else:
    raise NotImplementedError("Needs --100_times_dir option or --3_times_dirs, --50_times_dir, and --47_times_dir")

combined = defaultdict(set)

for count, td in enumerate(times_dirs):
    dg_dump_file = td + "/DG_state_plot.dump"
    if not OSU.check_file_exists(dg_dump_file):
        if OSU.check_file_exists(td + "results_except_draw.tgz"):
            print td + "results_except_draw.tgz: unpacking DG_state_plot.dump"
            OSU.system_command("tar -zxvf %sresults_except_draw.tgz ./DG_state_plot.dump" % (td))
            OSU.system_command("mv ./DG_state_plot.dump %s" % (td))
        else:
            raise IOError("results_except_draw.tgz not found in " + td)

    with open(dg_dump_file, "r") as f:
        print "Reading: " + dg_dump_file
        f.readline()  # throw away header
        for line in f:
            vars = line.split()
            str_key = "%s,%s" % (vars[0], vars[1])
            if vars[3] == "1" and vars[-1] == "1":
                combined[str_key].add(count)
            elif str_key not in combined:
                combined[str_key]

with open("%s/%scombined_DG_state_plot.dump" % (output_dir, file_prefix), "w") as f:
    f.write("nt\tDG\t%s\n" % ("\t".join([str(i) for i in range(1, 101)])))
    points = [k.split(",") for k in combined.keys()]
    points_sorted = sorted(points, key=lambda element: (float(element[0]), float(element[1])))
    for point in points_sorted:
        l = list(point) + ['0']*100
        for best_i in combined[",".join(point)]:
            l[best_i + 2] = '1'
        f.write("\t".join(l) + "\n")
