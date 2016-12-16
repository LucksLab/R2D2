import OSU
import VIU
from collections import defaultdict
import glob
import re
import math


def read_all_dbn_dirs(dbn_dirs):
    """
    Read dbns in multiple directories to pair correct lengths together.
    Returns a dictionary of .png files name that should have been made by R2D2 organized by length.
    """
    all_dbns = defaultdict(list)
    for dbn_dir in dbn_dirs:
        dbns = glob.glob(dbn_dir + "*.dbn")
        for dbn_f in dbns:
            with open(dbn_f, "r") as f:
                lines = f.readlines()
            length = len(lines[-1].split()[0])
            image = re.sub('(_mult\d+)?.dbn', '_structure.png', dbn_f)
            if image not in all_dbns[length] and OSU.check_file_exists(image):
                all_dbns[length].append(image)
    return all_dbns


if __name__ == "__main__":
    # read in arguments
    opts = OSU.getopts("", ["dbn_dirs=", "output_dir="])
    dbn_dirs = opts["--dbn_dirs"].split(",")
    width = 1200*len(dbn_dirs)
    output_dir = OSU.create_directory(opts["--output_dir"])

    # read dbns to pair correct lengths together
    all_dbns = read_all_dbn_dirs(dbn_dirs)

    # create images by horizontally concatenating previously made images from R2D2 output
    count = 0
    zero_padding = int(math.floor(math.log10(len(all_dbns))) + 1)
    for len in sorted(all_dbns):
        count += 1
        VIU.horizontal_image_concat("%s/%s.png"%(output_dir, str(count).zfill(zero_padding)), all_dbns[len])
        VIU.convert_center_resize("%s/%s.png"%(output_dir, str(count).zfill(zero_padding)), "%sx2800"%(width))
    print "ffmpeg -framerate 1 -i %s -vcodec mpeg4 -r 10 -s %s -pix_fmt yuv420p %s" % (output_dir + "/%%%dd.png" % (zero_padding), "%sx2800"%(width), output_dir + "/movie.mp4")
    OSU.system_command("ffmpeg -framerate 1 -i %s -vcodec mpeg4 -r 10 -s %s -pix_fmt yuv420p %s" % (output_dir + "/%%%dd.png" % (zero_padding), "%sx2800"%(width), output_dir + "/movie.mp4"))
