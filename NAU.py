"""
NAU (Nucleic Acid Utility) contains a variety of functions useful for manipulating RNA/DNA sequences and files

Author: This module was created by Kyle Watters in September, 2014, extended by Angela Yu 2014-2016.
Version: 0.0.0
"""

import sys
import getopt
import os.path
try:
    from collections import OrderedDict
except ImportError:
    from ordereddict import OrderedDict

name = "NAU.py"
operations = ['discover', 'rev', 'read']

help_message = '''

{0} (Nucleic Acid Utility) is a utility that contains a number of functions
(operations) that perform basic RNA/DNA manipulation

Usage:
    {0} [options] <operation> <operation input(s)>

Operations:   Choose an operation to perform

{1}      Determines the read length of a fastq file (1 input: file.fq or file.fastq)
{2}           Converts a DNA string to its reverse complement (1 input: DNA_string)
{3}          Converts a fasta-formatted file to a python dictionary of
              targets and sequences (1 input: file.fa or file.fasta)

Options:
-h,--help           brings up help for executing script
-f,--functions      brings up list of functions for import and use in other scripts
-v,--version        displays version number

'''.format(name, operations[0], operations[1], operations[2])


class Usage(Exception):
    def __init__(self,msg):
        self.msg = msg

def get_version():
    return "0.0.1"

class Params:
    def __init__(self):
        pass

    def parse_options(self, argv):

        try:
            opts, args = getopt.getopt(argv[1:],"hvq",["help","version","quiet"])

        except getopt.error, msg:
            raise Usage(msg)

        quiet = False

        for option, arg in opts:
            if option in ("-v", "--version"):
                print "%s v%s" % (name, get_version())
                exit(0)
            if option in ("-h", "--help"):
                raise Usage(help_message)
            if option in ("-q", "--quiet"):
                quiet = True

        if len(args) < 2 or (args[0] not in operations) \
                or (len(args) < 3 and args[0] == operations[3]):
            raise Usage(help_message)

        return args, quiet

    def check(self):
        pass


def read(fasta):
    #Check that the file is fasta
    extension = fasta.split(".")[-1]
    if extension != "fa" and extension != "fasta":
        print("File is not a fasta extension (fa,fasta)")

    #Check file exists
    elif not os.path.isfile(fasta):
        print("No such file as {0}".format(fasta))

    else:
        #Load the fasta file to parse
        targets = file(fasta, "r")

        lines = []
        for line in targets:
            read_line = line.strip()
            if read_line != "":
                if read_line[0] == ">":
                    read_line = read_line + " "
                lines.append(read_line)
        split_lines = ''.join(lines).split(">")[1:]
        all_targets = OrderedDict()
        for target in split_lines:
            #Takes the target name and sequence pairs, and splits them into a dictionary
            all_targets[target.split(" ")[0]] = target.split(" ")[1].strip()

        return all_targets


def rev(s):
    #This section was taken from Cole's code
    nuc_table = { 'A' : 'T',
                'T' : 'A',
                'C' : 'G',
                'G' : 'C',
                'U' : 'A',
                'a' : 't',
                't' : 'a',
                'c' : 'g',
                'g' : 'c',
                'u' : 'a',  }
    sl = list(s)

    try:
        rsl = [nuc_table[x] for x in sl]
    except KeyError, k:
        print >> sys.stderr, "Error: adapter sequences must contain only A,C,G,T,U"
        exit(1)
    rsl.reverse()

    return ''.join(rsl)


def discover(fastq):
    # Check that the file is fasta
    extension = fastq.split(".")[-1]
    if extension != "fq" and extension != "fastq":
        print("File is not a fasta extension (fq,fastq)")

    # Check file exists
    elif not os.path.isfile(fastq):
        print("No such file as {0}".format(fastq))

    else:
        # Detect the read length in the collapsed file
        fp = open(fastq)
        for i, line in enumerate(fp):
            if i == 1:
                read_len = len(line.strip())
                break
        fp.close()
        print("({0}) Read length is {1} nt.".format(name, read_len))

        return read_len


def make_seq(seq, seqfilename):
    seq = format_rna_string(seq)

    with open(seqfilename, "w") as new_seq:
        line_to_write = ';\n{0}\n%s1'.format(seqfilename) % (seq)
        new_seq.write(line_to_write)

    return seqfilename


def seq_to_dat_file(seqfile, datfile):
    with open(seqfile, 'r') as f:
        f.readline()
        name = f.readline()
        seq = f.readline()
        with open(datfile, 'w') as w:
            w.write("< " + name)
            w.write(seq[:-1])


def get_seq_from_ct(ctfile):
    with open(ctfile, 'r') as f:
        seqs = []
        seqcurr = []
        for line in f:
            vars = line.split()
            if len(vars) == 6 and "ENERGY" not in vars:
                seqcurr.append(format_rna_string(vars[1]))
            elif len(seqcurr) > 0:
                seqs.append("".join(seqcurr))
                seqcurr = []
        seqs.append("".join(seqcurr))  # last seq found not handled in loop
        return seqs


def format_rna_string(nt):
    """
    Capitalizes and replaces "T" with "U"
    """
    return nt.upper().replace('T', 'U')


def end_match_strip(seq, adapter):
    """
    Iteratively checks if end of seq matches the adapter and removes the match
    from the end of seq. If match is not found, the last character in adapter
    is removed until a match is found or all characters in adapter are removed.
    """
    for i in range(0, len(adapter)):
        substr = adapter[:-i]
        if i == 0:
            substr = adapter
        if seq.endswith(substr):
            return seq[:-len(substr)], len(adapter) - i
    return seq, -len(seq)


def main(argv=None,):

    params = Params()

    try:
        if argv is None:
            argv = sys.argv
            args, quiet = params.parse_options(argv)
            params.check()

        operation = args[0]
        input_arg = args[1]

        output = None
        if operation == operations[0]:
            output = discover(input_arg)
        elif operation == operations[1]:
            output = rev(input_arg)
        elif operation == operations[2]:
            output = read(input_arg)
        else:
            print("Operation not found")
            exit
        return output

    except Usage, err:
        print >> sys.stderr, sys.argv[0].split("/")[-1] + ": " + str(err.msg)
        print >> sys.stderr, ""
        return 2

if __name__ == "__main__":

    sys.exit(main())
