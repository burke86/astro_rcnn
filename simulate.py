import os
import time
import subprocess
from multiprocessing.dummy import Pool as ThreadPool

SET = 0

def bash(command,print_out=True):
    if print_out: print(command)
    try: return subprocess.call(command.split())
    except: return -1

def mask(line):
    set = SET
    OUT_DIR = os.path.abspath("../trainingset/")
    out_dir = os.path.abspath(os.path.join(OUT_DIR,"set_%d/" % set))
    if len(line.strip()) != 0:
        i = line.split()[1].split('.')[0]
        if line.split()[-3] == "star":
            obj_class = "star"
        else:
            obj_class = "gal"
        # Make a new catalog for each source in the image
        with open("./examples/obj%s" % i,"w+") as fi:
            fi.write("rightascension 0\ndeclination 0\nfilter 2\nvistime 90.0\nnsnap 1\nobshistid %s\n" % i)
            fi.write(line)
        # Now run PhoSim with this single object and no background to make mask
        bash("./phosim examples/obj%s -c examples/training_nobg -s 4S -i decam -t 1 -e 0" % i)
        out_to = os.path.abspath(os.path.join(out_dir,"%s%s.fits.gz" % (obj_class,i)))
        out_from = "./output/decam_e_%s_f2_4S_E000.fits.gz" % i
        while not os.path.exists(out_from):
            time.sleep(1)
        bash("mv %s %s" % (out_from, out_to))
        os.remove("./examples/obj%s" % i)


def simulate(set=0):
    os.chdir("./phosim_release")
    OUT_DIR = os.path.abspath("../trainingset/")
    # Run PhoSim on chip S4 with typical stars and galaxies
    # g
    bash("./phosim examples/maskrcnn_catalog_g -c examples/training -s 4S -i decam -t 128 -e 0")
    out_dir = os.path.abspath(os.path.join(OUT_DIR,"set_%d/" % set))
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
    out_to = os.path.abspath(os.path.join(out_dir,"img.fits.gz"))
    bash("mv ./output/decam_e_9999999_f2_4S_E000.fits.gz %s" % out_to)
    # r
    bash("./phosim examples/maskrcnn_catalog_r -c examples/training -s 4S -i decam -t 128 -e 0")
    out_dir = os.path.abspath(os.path.join(OUT_DIR,"set_%d/" % set))
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
    out_to = os.path.abspath(os.path.join(out_dir,"img.fits.gz"))
    bash("mv ./output/decam_e_9999999_f2_4S_E000.fits.gz %s" % out_to)
    # i
    bash("./phosim examples/maskrcnn_catalog_i -c examples/training -s 4S -i decam -t 128 -e 0")
    out_dir = os.path.abspath(os.path.join(OUT_DIR,"set_%d/" % set))
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
    out_to = os.path.abspath(os.path.join(out_dir,"img.fits.gz"))
    bash("mv ./output/decam_e_9999999_f2_4S_E000.fits.gz %s" % out_to)
    # Read the trimmed catalog for chip S4 and run on each object to generate mask
    with open("./work/trimcatalog_9999_4S.pars","r") as f:
        lines = []
        readlines = f.readlines()
        for i,line in enumerate(readlines):
            if len(line.strip()) != 0:
                obj_id = line.split()[1]
                if ".1" in obj_id: # 1st component of galaxy
                    # add second component
                    lines.append(line+"\n"+readlines[i+2])
                if ".2" in obj_id:
                    continue # already added 2nd component
                else: # star
                    lines.append(line)
        pool = ThreadPool(128)
        out = pool.map(mask, lines)
        pool.close()
        pool.join()
    # Remove everything in the work directory
    for f in os.listdir("./work"):
        os.remove(os.path.abspath(f))

if __name__ == "__main__":
    simulate(SET)
