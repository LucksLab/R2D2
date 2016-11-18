"""
Configuration for specific users during runtime.

config() only affects the current python process and changes will be lost once
the python process that calls this function exits. A more permanent way would
be to edit ~/.bash_profile (or ~/.bashrc) and define the environment variables.

Author: Angela M Yu, 2014-2016
Version: 0.0.1

Copyright (C) 2016  Julius B. Lucks, Angela M Yu.
All rights reserved.
Distributed under the terms of the GNU General Public License, see 'LICENSE'.
"""
import OSU


def config(user):
    """
    Configuration for each user
    """
    if user == "Julius":
        OSU.setenv("PATH", "/Users/jblucks/Software/RNAstructure/exe/")
        OSU.setenv("DATAPATH", "/Users/jblucks/Software/RNAstructure/data_tables/")
        OSU.setenv("CLASSPATH", "/Users/jblucks/Software/VARNA/VARNAv3-9.jar")
    elif user == "Angela":
        OSU.setenv("PATH", "/home/ayu/rotations/Lucks/RNAstructure/exe")
        OSU.setenv("DATAPATH", "/home/ayu/rotations/Lucks/RNAstructure/data_tables/")
        OSU.setenv("CLASSPATH", "/home/ayu/bin/VARNAv3-9.jar")
    elif user == "ICSE":
        OSU.setenv("PATH", ["/fs/europa/g_jbl/Software/Executables", "/usr/local/R/icse/bin", "/fs/home/amy35/GitHub/CoTrans_SHAPE-Seq"])
        OSU.setenv("DATAPATH", "/fs/europa/g_jbl/Software/RNAstructure/data_tables/")
        OSU.setenv("CLASSPATH", "/fs/europa/g_jbl/Software/Executables/VARNAv3-9.jar")
        OSU.setenv("LD_LIBRARY_PATH", "/fs/europa/g_jbl/Software/SHAPE-Seq/bin/lib:$LD_LIBRARY_PATH")
        OSU.setenv("BOOST_ROOT", "/fs/europa/g_jbl/Software/SHAPE-Seq/build/boost_1_49_0")
    elif user == "Quest":
        OSU.setenv("DATAPATH", "/projects/b1044/Software/src/RNAstructure/data_tables")
        OSU.setenv("CLASSPATH", "/projects/b1044/Software/bin/VARNAv3-9.jar")
        OSU.setenv("LD_LIBRARY_PATH", "/projects/b1044/Software/bin/SHAPE-Seq/bin/lib:$LD_LIBRARY_PATH")
        OSU.setenv("BOOST_ROOT", "/projects/b1044/Software/bin/SHAPE-Seq/build/boost_1_49_0")
    elif user == "Quest_R2D2":
        OSU.setenv("PATH", "/projects/b1044/Software/bin/GitHub/R2D2")
        OSU.setenv("DATAPATH", "/projects/b1044/Software/src/RNAstructure/data_tables")
        OSU.setenv("CLASSPATH", "/projects/b1044/Software/bin/VARNAv3-9.jar")
        OSU.setenv("LD_LIBRARY_PATH", "/projects/b1044/Software/bin/SHAPE-Seq/bin/lib:$LD_LIBRARY_PATH")
        OSU.setenv("BOOST_ROOT", "/projects/b1044/Software/bin/SHAPE-Seq/build/boost_1_49_0")

