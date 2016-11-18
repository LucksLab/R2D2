"""
Operating System Utilities (OSU)

Author: Angela M Yu, 2014-2016
Version: 0.0.1

Copyright (C) 2016  Julius B. Lucks, Angela M Yu.
All rights reserved.
Distributed under the terms of the GNU General Public License, see 'LICENSE'.
"""

import os
import glob
import sys
import getopt
from itertools import repeat, chain
import re
import cPickle as pickle


def setenv(env_var, paths_to_exe):
    """
    Sets $env_var environment variable, given a string or a list of strings.
    """
    oldpath = []
    if env_var in os.environ:
        oldpath = os.environ[env_var].split(":")
    if isinstance(paths_to_exe, list):
        os.environ[env_var] = ':'.join([p for p in paths_to_exe if p not in os.environ[env_var]])
    elif isinstance(paths_to_exe, str):
        if paths_to_exe not in oldpath:
            os.environ[env_var] = paths_to_exe
    else:
        raise TypeError("setenv input must be a string or list of strings")
    oldpath = ":".join(oldpath)
    if len(oldpath) > 0 and oldpath != os.environ[env_var]:
        os.environ[env_var] += ':' + oldpath


def create_directory(directory):
    """
        Create the output directory if it does not already exist.
        Return: The name of the directory.
    """
    if not os.path.exists(directory):
        os.makedirs(directory)
    return directory


def remove_files_with_ext(directory, extension):
    """ Remove all files in the specified directory with the extension """
    rmfiles = glob.glob(directory + "/*" + extension)
    for f in rmfiles:
        os.remove(f)


def remove_file(f):
    """ Removes a single file f """
    try:
        os.remove(f)
    except OSError:
        pass


def remove_files(files):
    """ Remove files in the list files """
    for f in files:
        remove_file(f)


def make_symbolic_link(source, link_name):
    """ Make a symbolic link of source to link_name """
    os.symlink(source, link_name)


def system_command(command):
    """ Execute a system command or list of commands """
    if isinstance(command, str):
        os.system(command)
    elif isinstance(command, list) and isinstance(command[0], str):
        for c in command:
            os.system(c)
    else:
        raise TypeError("OSU.system_command input must be a string or list of strings")


def getopts(short_arg_string, long_arg_list):
    """
    Returns a dictionary of command line arguments as defined by short_arg_string and long_arg_list
    """
    try:
        opts, args = getopt.getopt(sys.argv[1:], short_arg_string, long_arg_list)
    except getopt.GetoptError as err:
        print str(err)
        sys.exit(2)
    return dict(opts)


def flatten_list(list):
    """
    Flattens a list.
    From http://stackoverflow.com/questions/5286541/how-can-i-flatten-lists-without-splitting-strings
    """
    for x in list:
        if hasattr(x, '__iter__'):
            for y in flatten_list(x):
                yield y
        else:
            yield x


def load_all_pickles(indir):
    """
    Loads all pickles in a directory
    """
    pickles = {}
    infiles = glob.glob(indir + "/*.p")
    for f in infiles:
        fname = re.findall("[^/]*save_(.+).p$", f)
        if fname[0] != '':
            pickles[fname[0]] = pickle.load(open(f))
    return pickles


def increment_dict_counts(in_dict, update_dict):
    """
    Increments counts associated with the keys in in_dict (key -> int), based on update_dict.
    Assumes a collections.default_dict(int). Should consider extending to work with standard dictionaries as well as allow update_dict to be other data types (ex. list).
    """
    for k in update_dict:
        if k in in_dict:
            in_dict[k] += update_dict[k]
        else:
            in_dict[k] = update_dict[k]
    return in_dict


def ncycles(iterable, n):
    """
    Returns the sequence elements n times.
    From https://docs.python.org/2/library/itertools.html.
    """
    return chain.from_iterable(repeat(tuple(iterable), n))


def repeat_elements(list, i):
    """
    Returns a list containing each element of then list repeated i times.
    """
    return [x for x in list for j in range(i)]


def find_overlapping_string_indices(substr, string):
    """
    Searches string for all instance of a given substring and returns list of starting indices.
    Standard string and regex functions do not handle overlapping substring matches.
    """
    return [s.start() for s in re.finditer('(?=%s)' % (substr), string)]


def jsub_all_dir(dir):
    """
    Submits all *.sh files in directory to cluster
    """
    infiles = glob.glob(dir + "/*.sh")
    for f in infiles:
        system_command("/opt/voyager/nbs/bin/jsub " + f)


def check_file_exists(filename):
    """ Checks if file exists """
    return os.path.isfile(filename)


def get_dirname(filepath):
    """ Get directory from a file's path """
    return os.path.dirname(filepath)


def system_exit(comment=""):
    """ Exits with comment """
    sys.exit(comment)
