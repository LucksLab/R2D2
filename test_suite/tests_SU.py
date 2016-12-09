"""
Tests for SU.py in R2D2

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
import SU
<<<<<<< HEAD
=======
import LucksLabUtils_config
>>>>>>> JBL-code-review
from nose.tools import assert_equals
from numpy.testing import assert_array_almost_equal, assert_almost_equal, assert_equal
from numpy import array_equal
import os


def test_merge_labels():
    structs = [('1,1,1,1', "1"), ('1,1,1,1', "1"), ('1,1,1,1', "1_1"), ('2,2,2,2', "2"), ('3,3,3,3', "3,1_1")]
    merged_structs = SU.merge_labels(structs, False)
    print(str(merged_structs))
    assert_equals(set(merged_structs), set([('1,1,1,1', "1,1_1"), ('2,2,2,2', "2"), ('3,3,3,3', "3,1_1")]))


def test_merge_labels_toString():
    structs = [(['1','1','1','1'], "1"), (['1','1','1','1'], "1"), (['1','1','1','1'], "1_1"), (['2','2','2','2'], "2"), (['3','3','3','3'], "3,1_1")]
    merged_structs = SU.merge_labels(structs, True)
    print(str(merged_structs))
    assert_equals(set(merged_structs), set([('1,1,1,1', "1,1_1"), ('2,2,2,2', "2"), ('3,3,3,3', "3,1_1")]))


def test_calc_rho_from_theta_list():
    theta = [0.1, 0.2, 0.3, 0.4]
    assert_equals(SU.calc_rho_from_theta_list(theta), [0.4, 0.8, 1.2, 1.6])


def test_recalc_rhos():
    rhos = [0.4, 0.8, 1.2, 1.6]
    assert_array_almost_equal(SU.recalc_rhos(rhos, 1, 4), [0.8*3/3.6, 1.2*3/3.6, 1.6*3/3.6])
    assert_array_almost_equal(SU.recalc_rhos(rhos, 0, 3), [0.4*3/2.4, 0.8*3/2.4, 1.2*3/2.4])
    assert_array_almost_equal(SU.recalc_rhos(rhos, 0, 4), [0.4, 0.8, 1.2, 1.6])


def test_recalc_thetas():
    theta = [0.1, 0.2, 0.3, 0.4]
    assert_array_almost_equal(SU.recalc_thetas(theta, 0, 4), [0.1, 0.2, 0.3, 0.4])
    assert_array_almost_equal(SU.recalc_thetas(theta, 1, 4), [0.2/0.9, 0.3/0.9, 0.4/0.9])
    assert_array_almost_equal(SU.recalc_thetas(theta, 0, 3), [0.1/0.6, 0.2/0.6, 0.3/0.6])
    assert_array_almost_equal(SU.recalc_thetas(theta, 1, 3), [0.2/0.5, 0.3/0.5])


def test_invert_scale_rho_vec():
    rhos = [0.4, 0.8, 1.2, 1.6]
    assert_array_almost_equal(SU.invert_scale_rho_vec(rhos), [1-0.4/1.6, 1-0.8/1.6, 1-1.2/1.6, 0])
    assert_array_almost_equal(SU.invert_scale_rho_vec(rhos, 1), [0.6, 0.2, 0, 0])
    assert_array_almost_equal(SU.invert_scale_rho_vec(rhos, 1.5), [1-0.4/1.5, 1-0.8/1.5, 1-1.2/1.5, 0])


def test_scale_vec_avg1():
    rhos = [0.4, 0.8, 1.2, 2.6]
    assert_array_almost_equal(SU.scale_vec_avg1(rhos), [0.4*0.8, 0.8*0.8, 1.2*0.8, 2.6*0.8])
    assert_array_almost_equal(SU.scale_vec_avg1(rhos, 1), [0.4*4/3.2, 0.8*4/3.2, 1*4/3.2, 1*4/3.2])
    assert_array_almost_equal(SU.scale_vec_avg1([0]), [1])
    assert_array_almost_equal(SU.scale_vec_avg1([0], 0), [1])


def test_get_indices_rho_gt_c():
    rhos = [0.4, 0.8, 1.2, 2.6]
    assert_array_almost_equal(SU.get_indices_rho_gt_c(rhos, 2.6), [])
    assert_array_almost_equal(SU.get_indices_rho_gt_c(rhos, 1), [2, 3])
    assert_array_almost_equal(SU.get_indices_rho_gt_c(rhos, 2.6, True), [])
    assert_array_almost_equal(SU.get_indices_rho_gt_c(rhos, 1, True), [3, 4])


def test_ct_struct_to_binary_vec():
    ct = [5, 4, 0, 2, 1]
    ctl = [ct, [0, 3, 2]]
    assert_equals(SU.ct_struct_to_binary_vec(ct), [1, 1, 0, 1, 1])
    assert_equals(SU.ct_struct_to_binary_vec(ctl), [[1, 1, 0, 1, 1], [0, 1, 1]])


def test_ct_struct_to_binary_mat():
    ct = [2, 1, 0]
    ctl = [ct, [0, 3, 2]]
    assert_equals(SU.ct_struct_to_binary_mat(ct), [[0,1,0], [1,0,0], [0,0,0]])
    assert_equals(SU.ct_struct_to_binary_mat(ctl), [[[0,1,0], [1,0,0], [0,0,0]], [[0,0,0], [0,0,1], [0,1,0]]])


def test_cap_rho_or_ct_list():
    ct = [2, 1, 0]
    assert_equals(SU.cap_rho_or_ct_list(ct), [2,1,0])
    assert_equals(SU.cap_rho_or_ct_list(ct, 1), [1,1,0])


def test_calc_bp_distance_vector_weighted():
    struct = [1, 1, 0]
    rho = [1, 2, 0]
    assert_equals(SU.calc_bp_distance_vector_weighted(struct, rho), 0.5)
    assert_equals(SU.calc_bp_distance_vector_weighted(struct, rho, scaling_func="D"), 0.5)
    assert_equals(SU.calc_bp_distance_vector_weighted(struct, rho, scaling_func="K"), 0.5)
    assert_equals(SU.calc_bp_distance_vector_weighted(struct, rho, scaling_func="U"), 0.5)
    assert_equals(SU.calc_bp_distance_vector_weighted(struct, rho, scaling_func="D", max_r=0.9), 0.5*(0.1+1.1))
    assert_equals(SU.calc_bp_distance_vector_weighted(struct, rho, scaling_func="D", invert_struct=True), 3*0.5+0.5)
    assert_equals(SU.calc_bp_distance_vector_weighted(struct, rho, scaling_func="D", paired_weight=1), 1)

    struct = [1, 0, 1] # paired, unpaired, paired
    rho = SU.invert_scale_rho_vec([2, 0, 0], 1.8)

    assert_equals(SU.calc_bp_distance_vector_weighted(struct, rho, scaling_func="D", paired_weight=0.6), 1)
    assert_equals(SU.calc_bp_distance_vector_weighted(struct, rho, scaling_func="D", paired_weight=0.9), 1)
    assert_equals(SU.calc_bp_distance_vector_weighted(struct, rho, scaling_func="D", paired_weight=0.1), 1)
    assert_equals(SU.calc_bp_distance_vector_weighted(struct, rho, scaling_func="D", paired_weight=0.3), 1)

    struct = [1, 0, 1] # paired, unpaired, paired
    rho = SU.invert_scale_rho_vec([2, 0.1, 0], 1.8)

    assert_almost_equal(SU.calc_bp_distance_vector_weighted(struct, rho, scaling_func="D", paired_weight=0.6), 0.6+0.4*(1-0.1/1.8))
    assert_almost_equal(SU.calc_bp_distance_vector_weighted(struct, rho, scaling_func="D", paired_weight=0.9), 0.9+0.1*(1-0.1/1.8))
    assert_almost_equal(SU.calc_bp_distance_vector_weighted(struct, rho, scaling_func="D", paired_weight=0.1), 0.1+0.9*(1-0.1/1.8))
    assert_almost_equal(SU.calc_bp_distance_vector_weighted(struct, rho, scaling_func="D", paired_weight=0.3), 0.3+0.7*(1-0.1/1.8))


def test_calc_bp_distance_matrix():
    react_mat = [[0,1,0], [0,0,0], [0,0,0]]
    ct_mat = [[0,0,0], [0,0,0], [0,0,0]]
    assert_equals(SU.calc_bp_distance_matrix(react_mat, ct_mat), 1)
    assert_equals(SU.calc_bp_distance_matrix(react_mat, ct_mat, endoff=1), 0)

    react_mat = [[0,1,0], [1,0,0], [0,0,0]]
    ct_mat = [[0,0,0], [0,0,0], [0,0,0]]
    assert_equals(SU.calc_bp_distance_matrix(react_mat, ct_mat), 1)
    assert_equals(SU.calc_bp_distance_matrix(react_mat, ct_mat, endoff=1), 0)


def test_calc_distances_bt_matrices():
    react_mat = SU.ct_struct_to_binary_mat([2,1,0])
    ct_mat = SU.ct_struct_to_binary_mat([0,0,0])
    assert array_equal(SU.calc_distances_bt_matrices([react_mat, ct_mat]), [[0, 1], [1, 0]]), "test_calc_distances_bt_matrices: Failed"

    react_mat = SU.ct_struct_to_binary_mat([0,0,2,1])
    ct_mat = SU.ct_struct_to_binary_mat([0,0,0,0])
    assert array_equal(SU.calc_distances_bt_matrices([react_mat, ct_mat], 2), [[0, 0], [0, 0]]), "test_calc_distances_bt_matrices: Failed endoff=2"


def test_calc_benchmark_statistics_matrix():
    react_mat = SU.ct_struct_to_binary_mat([0,0,0])
    ct_mat = SU.ct_struct_to_binary_mat([2,1,0])
    assert_equal(SU.calc_benchmark_statistics_matrix(react_mat, ct_mat), {"Sensitivity":0., "PPV":float('nan'), "F_score":0.})

    react_mat = SU.ct_struct_to_binary_mat([0,0,0])
    ct_mat = SU.ct_struct_to_binary_mat([0,0,0])
    assert_equal(SU.calc_benchmark_statistics_matrix(react_mat, ct_mat), {"F_score":float('nan'), "Sensitivity":float('nan'), "PPV":float('nan')})

    react_mat = SU.ct_struct_to_binary_mat([2,1,0])
    ct_mat = SU.ct_struct_to_binary_mat([2,1,0])
    assert_equal(SU.calc_benchmark_statistics_matrix(react_mat, ct_mat), {"F_score":1., "Sensitivity":1., "PPV":1.})

    react_mat = SU.ct_struct_to_binary_mat([2,1,4,3])
    ct_mat = SU.ct_struct_to_binary_mat([2,1,0,0])
    assert_equal(SU.calc_benchmark_statistics_matrix(react_mat, ct_mat), {"F_score":2/3., "Sensitivity":1., "PPV":0.5})


def remove_file(f):
    try:
        os.remove(f)
    except OSError:
        pass

class TestSU:
    @classmethod
    def setup_class(cls):
        cls.input_dir = "."
        cls.outputfile = "test"
        cls.reactivity_file = "test_reactivities.txt"
        cls.adapterseq = "CTGACTCGGGCACCAAGG"
        cls.theta_file = "test_theta.theta"
        cls.ctfile_1 = "test_1.ct"
        cls.ctfile_2 = "test_2.ct"


    @classmethod
    def teardown_class(cls):
        remove_file("test")
        remove_file("test.rho")
        remove_file("test.seq")
        remove_file("test.theta")


    def test_parse_reactivity_rho(self):
        pos, rho, theta_cut, rho_cut, seq_cut, rc_flag, rc_sum, untreated_sum, treated_sum = SU.parse_reactivity_rho(self.reactivity_file, self.adapterseq, self.outputfile)
        assert_equals(pos, [str(a) for a in range(1,21)])
        assert_array_almost_equal(rho, [1]*20)
        assert_equals(theta_cut, ["0.05"]*20)
        assert_array_almost_equal(rho_cut, [1.0]*20)
        assert_equals(seq_cut, "AUCGGGGGCUCUGUUGGUUC")
        assert_equals(rc_flag, 1)
        assert_equals(rc_sum, 2050+2*19+19)
        assert_equals(untreated_sum, 1010+19)
        assert_equals(treated_sum, 1040+2*19)

        pos, rho, theta_cut, rho_cut, seq_cut, rc_flag, rc_sum, untreated_sum, treated_sum = SU.parse_reactivity_rho(self.reactivity_file, self.adapterseq, self.outputfile, endcut=-14)
        assert_equals(pos, [str(a) for a in range(1,7)])
        assert_array_almost_equal(rho, [1]*20)
        assert_equals(theta_cut, [str(1/6.0)]*6)
        assert_array_almost_equal(rho_cut, [1.0]*6)
        assert_equals(seq_cut, "AUCGGG")
        assert_equals(rc_flag, 1)
        assert_equals(rc_sum, 2050+2*5+5)
        assert_equals(untreated_sum, 1010+5)
        assert_equals(treated_sum, 1040+2*5)


    #def test_reactivities_to_rho_file(self):

    #reactivities_to_reads_files


    def test_calc_rho_from_theta(self):
        assert_equals(SU.calc_rho_from_theta(self.theta_file, "test"), ['1.0', '1.0', '2.0', '0.0'])


    #rhos_list_to_file

    #ct_list_to_file


    def test_cts_to_file(self):
       cts = SU.get_ct_structs(self.ctfile_1)
       SU.cts_to_file(cts, "A", self.outputfile)
       with open(self.ctfile_1, "r") as f:
           lines_orig = f.readlines()
       with open(self.outputfile, "r") as f:
           lines = f.readlines()
       for l in zip(lines_orig, lines):
           orig = l[0].split()
           if "/" in orig[1]:
               assert_equals(orig[0], l[1].split()[0])
           else:
               assert_equals(orig, l[1].split())

       cts = SU.get_ct_structs(self.ctfile_2)
       SU.cts_to_file(cts, "AUCGGGGG", self.outputfile)
       with open(self.ctfile_2, "r") as f:
           lines_orig = f.readlines()
       with open(self.outputfile, "r") as f:
           lines = f.readlines()
       for l in zip(lines_orig, lines):
           orig = l[0].split()
           if "/" in orig[1]:
               assert_equals(orig[0], l[1].split()[0])
           else:
               assert_equals(orig, l[1].split())


    #make_constraint_file


    def test_get_free_energy_efn2(self):
        LucksLabUtils_config.config("Quest_R2D2")
        SU.runRNAstructure_efn2(self.ctfile_1, self.outputfile)
        energy = SU.get_free_energy_efn2(self.outputfile)
        assert_equals([0], energy)

        SU.runRNAstructure_efn2(self.ctfile_2, self.outputfile)
        energy = SU.get_free_energy_efn2(self.outputfile)
        assert_equals([5.7, 0.0, 4.1, 3.9], energy)


    def test_get_ct_structs(self):
        cts = SU.get_ct_structs(self.ctfile_1)
        assert_equals([['0']], cts)

        cts = SU.get_ct_structs(self.ctfile_2)
        zeros = ["0"] * 8
        orig = [zeros[:], zeros[:], zeros[:], zeros[:]]
        orig[0][1] = '8'
        orig[0][-1] = '2'
        orig[2][1] = '8'
        orig[2][-1] = '2'
        orig[2][2] = '7'
        orig[2][-2] = '3'
        orig[3][2] = '7'
        orig[3][-2] = '3'
        assert_equals(orig, cts)


    #ct_file_to_struct_file
