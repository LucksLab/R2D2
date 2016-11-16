import SU
import OSU
import PAU
import NAU
import LucksLabUtils_config
import glob

# setup environment variables specific to the ICSE cluster at Cornell
LucksLabUtils_config.config("Quest_R2D2")
opts = OSU.getopts("o:c:r:p:", ["shape_intercept=", "shape_slope="])
print opts

reactivities_files = glob.glob(opts['-r'])
crystal_files = glob.glob(opts['-c'])
output_dir = OSU.create_directory(opts['-o'])
num_proc = int(opts['-p']) if '-p' in opts else 1
shape_intercept = float(opts["--shape_intercept"]) \
                    if "--shape_intercept" in opts else -0.3
shape_slope = float(opts["--shape_slope"]) if "--shape_slope" in opts else 1.1

reactivities = PAU.parse_input_panels(reactivities_files, output_dir)
crystals = PAU.parse_input_panels(crystal_files, output_dir)

if set(crystals.keys()) != set([rk.split('-')[0] for rk in reactivities.keys()]):
    raise Exception("Keys of crystal structures and reactivities not equal")

print crystals.keys()
print reactivities.keys()

for k, v in reactivities.iteritems():
    react_rhos = v[0]
    react_seq = v[1]
    ck = k.split('-')[0]  # corresponding crystal key from the reactivity key
    cryst_seq = crystals[ck][1]
    ind_of_match = react_seq.index(cryst_seq)
    end_of_match = ind_of_match + len(cryst_seq)
    # renormalize rhos based on matching regions
    react_rhos = SU.recalc_rhos(react_rhos, ind_of_match, end_of_match)

    reactivities[k][0] = react_rhos
    reactivities[k][1] = react_seq[ind_of_match:end_of_match]
    SU.rhos_list_to_file(react_rhos, reactivities[k][2]+".rho")
    seqfile = NAU.make_seq(reactivities[k][1], reactivities[k][2]+".seq")
    SU.runRNAstructure_fold(seqfile, reactivities[k][2]+".ct", shapefile=reactivities[k][2]+".rho", p=num_proc, shape_intercept=shape_intercept, shape_slope=shape_slope)
    SU.runRNAstructure_CircleCompare(reactivities[k][2]+".ct", crystals[ck][3], reactivities[k][2]+".ps")
    OSU.system_command("convert %s.ps %s.jpg" % (reactivities[k][2], reactivities[k][2]))
    with open(reactivities[k][2]+".stats", "w") as f:
        react_mat = SU.ct_struct_to_binary_mat(SU.get_ct_structs(reactivities[k][2]+".ct"))[0]
        f.write(str(SU.calc_benchmark_statistics_matrix(react_mat, crystals[ck][2])))
