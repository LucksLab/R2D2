"""
Video and Image Utilities (VIU.py)

Version: 0.0.0
Author: Angela M Yu, 2014-2016
"""

import OSU
import os
import NAU
import SU
import math
import glob


def vertical_image_concat(outfile, images):
    """
    Vertical concatenate images to a new image file.
    """
    OSU.system_command("convert -append " + " ".join(images) + " " + outfile)


def convert_center_resize(image, res):
    """
    Centers image and resizes.
    """
    OSU.system_command("convert %s -background none -gravity Center -resize %s %s.temp" % (image, res, image))
    os.rename(image + ".temp", image)


def convert_center(image):
    """
    Centers image.
    """
    OSU.system_command("convert %s -background none -gravity Center %s.temp" % (image, image))
    os.rename(image + ".temp", image)


def generate_movie(regex, outfile):
    """
    Generate a movie with images as described by regex.
    """
    try:
        print "ffmpeg -r 1 -i " + regex + " -vcodec mpeg4 -b 800k -r 10 -s 1200x2800 -pix_fmt yuv420p " + outfile
        OSU.system_command("ffmpeg -r 1 -i " + regex + " -vcodec mpeg4 -b 800k -r 10 -s 1200x2800 -pix_fmt yuv420p " + outfile)
    except:
        print "ffmpeg -framerate 1 -i " + regex + " -c:v libx264 -r 10 -s 1200x2800 -pix_fmt yuv420p " + outfile
        OSU.system_command("ffmpeg -framerate 1 -i " + regex + " -c:v libx264 -r 10 -s 1200x2800 -pix_fmt yuv420p " + outfile)


def generate_MFE_CoTrans_movie(seq, seq_start, seq_end, outdir, thetasdir=""):
    """
    Generate co-transcriptional MFE folding movie.
    """
    OSU.create_directory(outdir)
    OSU.create_directory(outdir + "/seq/")
    OSU.create_directory(outdir + "/ct/")
    zero_padding = int(math.floor(math.log10(seq_end)) + 1)
    varna_num = 0
    rhos = {}
    if thetasdir != "":
        # should use theta files to generate rhos because the adapter is already trimmed off
        for rf in glob.glob(thetasdir+"/*.theta"):
            rho = SU.calc_rho_from_theta(rf, "temp.rho")
            rhos[len(rho)] = rho  # add in rho file here

        OSU.remove_file("temp.rho")
    with open("VIU.test.txt", 'w') as vtf:
        vtf.write(str(rhos))
    for seqi in range(seq_start+1, seq_end+1):
        varna_num += 1
        if seqi in rhos:
            rho_varna = "\"" + ";".join(rhos[seqi]+(["-1"]*(seq_end-seqi))) + "\""
        else:
            rho_varna = "\"" + ";".join(["-1"]*(seq_end)) + "\""
        seqf = outdir + "/seq/" + str(seqi) + ".seq"
        ctf = outdir + "/ct/" + str(seqi) + ".ct"
        NAU.make_seq(seq[seq_start:seqi], seqf)
        SU.runRNAstructure_fold(seqf, ctf)
        SU.run_ct2dot(ctf, 0, "temp.dbn")
        OSU.system_command("sed '$s/$/&%s/' temp.dbn > temp_ext.dbn " % ("." * (seq_end - seqi)))
        run_VARNA("temp_ext.dbn", outdir + str(varna_num).zfill(zero_padding) + "_structure.png", rho_varna)
        convert_center_resize(outdir + str(varna_num).zfill(zero_padding) + "_structure.png", "1440x2500")
    OSU.remove_file(outdir + "temp.dbn")
    OSU.remove_file(outdir + "temp_ext.dbn")
    OSU.system_command("cp %s %s" % (outdir + str(varna_num).zfill(zero_padding) + "_structure.png", outdir + str(varna_num+1).zfill(zero_padding) + "_structure.png"))
    generate_movie(outdir + "%%%dd_structure.png" % (zero_padding), outdir + "/movie.mp4")


def run_VARNA(dbnfile, outfile, SHAPE_vals):
    """
    Run VARNA to generate .jpg file of secondary structure from .dbn file
    and SHAPE reactivities.
    """
    mapcolors = '"0:#FFFFFF;0.75:#FFCC33;1.3:#FF9933;2.0:#FF6633;5.0:#FF0000"'
    command = 'java fr.orsay.lri.varna.applications.VARNAcmd -i %s -o %s -resolution "3.0" -title "" -titleSize 10  -colorMap %s -colorMapStyle %s  -colorMapMin "0.0" -colorMapMax "5.0"' % (dbnfile, outfile, SHAPE_vals, mapcolors)
    #command += ' > /dev/null 2>&1'
    print command
    os.system(command)
