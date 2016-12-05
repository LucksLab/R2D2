"""
Tests for NAU.py in R2D2

Usage: Runs with 'nosetests' when called within the 'test_suite' directory
Ex.
[ayd0188@quser11 test_suite]$ nosetests
.........................
----------------------------------------------------------------------
Ran 25 tests in 0.110s

OK


Version: 0.0.1
Author: Angela M Yu, 2016

Copyright (C) 2016  Julius B. Lucks, Angela M Yu.
All rights reserved.
Distributed under the terms of the GNU General Public License, see 'LICENSE'.
"""

import sys
sys.path.insert(0,'..')
import NAU
from nose.tools import assert_equals
import os


def test_rev():
    seq = "AUCGaucg"
    assert_equals(NAU.rev(seq), "cgatCGAT")


def test_format_rna_string():
    seq = "AUCGaucgATCGatcg"
    assert_equals(NAU.format_rna_string(seq), "AUCGAUCGAUCGAUCG")


def test_end_match_strip():
    adapter = "AUCG"
    seq = "AUCGAUCG"
    assert_equals(NAU.end_match_strip(seq, adapter), ("AUCG", 4))

    adapter = "UUCG"
    seq = "AAAA"
    assert_equals(NAU.end_match_strip(seq, adapter), ("AAAA", -4))

    adapter = ""
    seq = "AUCGAUCG"
    assert_equals(NAU.end_match_strip(seq, adapter), ("AUCGAUCG", -8))

    adapter = "AUCG"
    seq = ""
    assert_equals(NAU.end_match_strip(seq, adapter), ("", 0))

    adapter = ""
    seq = ""
    assert_equals(NAU.end_match_strip(seq, adapter), ("", 0))


def remove_file(f):
    try:
        os.remove(f)
    except OSError:
        pass


class TestNAU:
    @classmethod
    def setup_class(cls):
        cls.seq = "CTGACTCGGGCACCAAGG"
        cls.seqfilename = "test.seq"
        cls.datfile = "test.dat"
        cls.ctfile_1 = "test_1.ct"
        cls.ctfile_2 = "test_2.ct"
        cls.fasta_1 = "test_1.fa"
        cls.fasta_2 = "test_2.fasta"
        cls.fastq_1 = "test_1.fq"
        cls.fastq_2 = "test_2.fastq"

    @classmethod
    def teardown_class(cls):
        remove_file(cls.seqfilename)
        remove_file(cls.datfile)


    def test_make_seq(self):
        seqfilename = NAU.make_seq(self.seq, self.seqfilename)
        assert_equals(seqfilename, self.seqfilename)
        assert_equals(os.path.isfile(seqfilename), True)
        with open(seqfilename, "r") as f:
            lines = f.readlines()
        assert_equals(lines[0], ";\n")
        assert_equals(lines[1], seqfilename + "\n")
        assert_equals(lines[2], NAU.format_rna_string(self.seq) + "1")


    def test_seq_to_dat_file(self):
        seqfile = seqfilename = NAU.make_seq(self.seq, self.seqfilename)
        NAU.seq_to_dat_file(self.seqfilename, self.datfile)
        assert_equals(os.path.isfile(self.datfile), True)
        with open(self.datfile, "r") as f:
            lines = f.readlines()
        assert_equals(lines[0], "< " + seqfilename + "\n")
        assert_equals(lines[1], NAU.format_rna_string(self.seq))


    def test_get_seq_from_ct(self):
        seqs = NAU.get_seq_from_ct(self.ctfile_1)
        assert_equals(seqs[0], "A")

        seqs = NAU.get_seq_from_ct(self.ctfile_2)
        assert_equals(seqs, ["AUCGGGGG"]*4)


    def test_read(self):
        assert_equals(NAU.read("test.no"), None)
        assert_equals(NAU.read("does_not_exist.fa"), None)
        assert_equals(NAU.read(self.fasta_1), {"1": "AUCG"})
        assert_equals(NAU.read(self.fasta_2), {"1": "AUCGA", "2": "GUCGG"})


    def test_discover(self):
        assert_equals(NAU.discover("test.no"), None)
        assert_equals(NAU.discover("does_not_exist.fq"), None)
        assert_equals(NAU.discover(self.fastq_1), 4)
        assert_equals(NAU.discover(self.fastq_2), 5)
