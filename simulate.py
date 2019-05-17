import os
import subprocess

NUM_IMAGES = 1

def bash(command,print_out=True):
    if print_out: print(command)
    try: return subprocess.call(command.split())
    except: return -1

def simulate(set=0):
    os.chdir("./phosim_release")
    # Run PhoSim on chip S4 with typical stars and galaxies
    bash("./phosim examples/maskrcnn_catalog -c examples/training_nobg -s 4S -i decam")
    if not os.path.exists("./output/set_%d" % set):
        os.mkdir("./output/set_%d" % set)
    bash("mv ./output/decam_e_9999_f2_4S_E000.fits.gz ./output/set_%d/img.fits.gz" % set)
    # Read the trimmed catalog for chip S4
    with open("./work/trimcatalog_9999_4S.pars","r") as f:
        for line in f.readlines():
            if len(line.strip()) != 0:
                i = line.split()[1]
                obj_class = line.split()[-3]
                # Make a new catalog for each source in the image
                with open("./examples/obj%s" % i,"w+") as fi:
                    fi.write("rightascension 0\ndeclination 0\nfilter 2\nvistime 90.0\nnsnap 1\n")
                    fi.write(line)
                # Now run PhoSim with this single object and no background to make mask
                bash("./phosim examples/obj%s -c examples/training_nobg -s 4S -i decam" % i)
                os.remove("./examples/obj%s" % i)
                if not os.path.exists("./output/set_%d" % set):
                    os.mkdir("./output/set_%d" % set)
                bash("mv ./output/decam_e_9999_f2_4S_E000.fits.gz ./output/set_%d/%s%s.fits.gz" % (set,obj_class,i))
    # Remove everything in the work directory
    for f in os.listdir("./phosim_release/work"):
        bash("rm %s" % os.abspath(f))


if __name__ == "__main__":
    for i in range(NUM_IMAGES):
        simulate(i)
