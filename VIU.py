"""
Video and Image Utilities (VIU.py)

Version: 0.0.1
Author: Angela M Yu, 2014-2016

Copyright (C) 2016  Julius B. Lucks, Angela M Yu.
All rights reserved.
Distributed under the terms of the GNU General Public License, see 'LICENSE'.
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


def horizontal_image_concat(outfile, images):
    """
    Horizontal concatenate images to a new image file.
    """
    OSU.system_command("convert +append " + " ".join(images) + " " + outfile)


def convert_center_resize(image, res):
    """
    Centers image and resizes.
    """
    try:
        print "convert %s -alpha discrete -blur 0x1 -background none -gravity Center -extent %s %s.temp" % (image, res, image)
        OSU.system_command("convert %s -alpha discrete -blur 0x1 -background none -gravity Center -extent %s %s.temp" % (image, res, image))
        os.rename(image + ".temp", image)
    except:
        print "convert %s -background none -gravity Center -extent %s %s.temp" % (image, res, image)
        OSU.system_command("convert %s -background none -gravity Center -extent %s %s.temp" % (image, res, image))
        os.rename(image + ".temp", image)


def convert_center(image):
    """
    Centers image.
    """
    OSU.system_command("convert %s -background none -gravity Center %s.temp" % (image, image))
    os.rename(image + ".temp", image)


def generate_movie(regex, outfile, size="1200x2800"):
    """
    Generate a movie with images as described by regex.
    """
    if size != "":
        try:
            print "ffmpeg -r 1 -i " + regex + " -vcodec mpeg4 -b 800k -r 10 -s " + size + " -pix_fmt yuv420p " + outfile
            OSU.system_command("ffmpeg -r 1 -i " + regex + " -vcodec mpeg4 -b 800k -r 10 -s " + size + " -pix_fmt yuv420p " + outfile)
        except:
            print "ffmpeg -framerate 1 -i " + regex + " -c:v libx264 -r 10 -s " + size + " -pix_fmt yuv420p " + outfile
            OSU.system_command("ffmpeg -framerate 1 -i " + regex + " -c:v libx264 -r 10 -s " + size + " -pix_fmt yuv420p " + outfile)
    else:
        print "ffmpeg -framerate 1 -i " + regex + " -vcodec mpeg4 -b 800k -r 10 -pix_fmt yuv420p " + outfile
        OSU.system_command("ffmpeg -framerate 1 -i " + regex + " -vcodec mpeg4 -b 800k -r 10 -pix_fmt yuv420p " + outfile)


def generate_MFE_CoTrans_movie(seq, outdir, seq_start=-1, seq_end=-1, rhos_dir="", SHAPE_direct=False):
    """
    Generate co-transcriptional MFE folding movie.
    Options to start and end at specific lengths, seq_start and seq_end respectively.
    Can overlay rho reactivities if given a directory with .rho files corresponding to the sequence.
    """
    OSU.create_directory(outdir)
    OSU.create_directory(outdir + "/seq/")
    OSU.create_directory(outdir + "/ct/")
    if seq_start == -1:
        seq_start = 0
    if seq_end == -1:
        seq_end = len(seq)
    else:
        seq_end += 1
    zero_padding = int(math.floor(math.log10(seq_end)) + 1)
    varna_num = 0
    rhos = {}
    if rhos_dir != "":
        # reads through .rho files found in rhos_dir
        for rf in glob.glob(rhos_dir+"/*.rho"):
            # read in each rho reactivitiy spectra
            with open(rf, "r") as f:
                rho = [line.split()[1] for line in f.readlines()]
                rhos[len(rho)] = [rho, rf]  # add in rho file here

    for seqi in range(seq_start+1, seq_end+1):
        if seqi in rhos:
            rho_varna = "\"" + ";".join(rhos[seqi][0]+(["-1"]*(seq_end-seqi))) + "\""
        else:
            rho_varna = "\"" + ";".join(["-1"]*(seq_end)) + "\""
        seqf = outdir + "/seq/" + str(seqi) + ".seq"
        ctf = outdir + "/ct/" + str(seqi) + ".ct"
        NAU.make_seq(seq[seq_start:seqi], seqf)
        if SHAPE_direct and seqi in rhos:
            SU.runRNAstructure_fold(seqf, ctf, rhos[seqi][1])
        elif SHAPE_direct:
            continue
        else:
            SU.runRNAstructure_fold(seqf, ctf)
        SU.run_ct2dot(ctf, 0, "temp.dbn")
        OSU.system_command("sed '$s/$/&%s/' temp.dbn > temp_ext.dbn " % ("." * (seq_end - seqi)))
        varna_num += 1
        run_VARNA("temp_ext.dbn", outdir + str(varna_num).zfill(zero_padding) + "_structure.png", rho_varna)
        convert_center_resize(outdir + str(varna_num).zfill(zero_padding) + "_structure.png", "1440x2000")
    OSU.remove_file(outdir + "temp.dbn")
    OSU.remove_file(outdir + "temp_ext.dbn")
    generate_movie(outdir + "%%%dd_structure.png" % (zero_padding), outdir + "/movie.mp4", "")


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
