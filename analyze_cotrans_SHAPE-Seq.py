"""
analyze_cotrans_SHAPE-Seq.py calls R2D2
See examples/run_CoTrans_example.sh for an example of how to use this code.

Author: Angela M Yu, 2014-2016
Version: 0.0.1

Copyright (C) 2016  Julius B. Lucks and Angela M Yu.
All rights reserved.
Distributed under the terms of the GNU General Public License, see 'LICENSE'.
"""

import R2D2
import LucksLabUtils_config
import OSU

import ipdb; ipdb.set_trace() #JBL- entering debugging here
LucksLabUtils_config.config("Quest_R2D2")
OSU.system_command("echo $PATH")
OSU.system_command("echo $CLASSPATH")

opts = OSU.getopts("", ["in_dir=", "out_dir=", "adapter=", "p=", "e=", "endcut=", "constrained_c=", "scale_rho_max=", "draw_all=", "most_count_tie_break=", "scaling_func=", "weight_paired=", "cap_rhos=", "pol_fp="])
print opts

# This specifically calls R2D2.R2D2() assuming the user has specified the arguments:
# in_dir, out_dir, adapter, e, endcut, constrained_c, scale_rho_max, draw_all, most_count_tie_break, scaling_func, weight_paired, cap_rhos, pol_fp
# Only in_dir, out_dir, and adapter are truly required to run R2D2.R2D2(). Default values for the other parameters are set within R2D2.py.

cotrans = R2D2.R2D2(opts['--in_dir'], opts['--out_dir'], opts['--adapter'], p=int(opts['--p']), e=int(opts['--e']), endcut=int(opts['--endcut']), constrained_c=float(opts['--constrained_c']), scale_rho_max=float(opts['--scale_rho_max']), draw_all=bool(opts['--draw_all'] == "True"), most_count_tie_break=bool(opts['--most_count_tie_break'] == "True"), scaling_func=opts['--scaling_func'], weight_paired=float(opts['--weight_paired']), cap_rhos=bool(opts["--cap_rhos"]=="True"), pol_fp=int(opts['--pol_fp']))
