import os
import time
import subprocess
from multiprocessing.dummy import Pool as ThreadPool

os.chdir("./phosim_release")

def bash(command,print_out=True):
    if print_out: print(command)
    try: return subprocess.call(command.split())
    except: return -1

class PhoSimSet:

    def __init__(self,set_num):
        self.set = set_num

    def mask(self,line):
        OUT_DIR = os.path.abspath("../trainingset/")
        out_dir = os.path.abspath(os.path.join(OUT_DIR,"set_%d/" % self.set))
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
                fi.write("\n")
            print(i)
            print(line)
            # Now run PhoSim with this single object and no background to make mask
            """bash("./phosim examples/obj%s -c examples/training_nobg -i decam -t 1 -e 0" % i)
            out_to = os.path.abspath(os.path.join(out_dir,"%s%s.fits.gz" % (obj_class,i)))
            out_from = "./output/decam_e_%s_f2_4S_E000.fits.gz" % i
            wait = 0
            while not os.path.exists(out_from):
                time.sleep(2)
                wait += 1
                print("wait")
                if wait == 10: break
            try:
                bash("mv %s %s" % (out_from, out_to))
                os.remove("./examples/obj%s" % i)
            except: # off chip
            """
        return

    def img(self,band):
        OUT_DIR = os.path.abspath("../trainingset/")
        out_dir = os.path.abspath(os.path.join(OUT_DIR,"set_%d/" % self.set))
        out_to = os.path.abspath(os.path.join(out_dir,"img_%s.fits.gz" % band))
        bash("./phosim examples/maskrcnn_catalog_%s -c examples/training -i decam -t 42 -e 0" % band)
        if band == 'g':
            bash("mv ./output/decam_e_9999999_f1_4S_E000.fits.gz %s" % out_to)
        elif band == 'r':
            bash("mv ./output/decam_e_9999998_f2_4S_E000.fits.gz %s" % out_to)
        elif band == 'i':
            bash("mv ./output/decam_e_9999997_f3_4S_E000.fits.gz %s" % out_to)
        return

    def simulate(self):
        #os.chdir("./phosim_release")
        OUT_DIR = os.path.abspath("../trainingset/")
        out_dir = os.path.abspath(os.path.join(OUT_DIR,"set_%d/" % self.set))
        if os.path.exists(out_dir):
            print("Error: Set %d already exists." % self.set)
            return
        os.mkdir(out_dir)
        # Run PhoSim on chip S4 with typical stars and galaxies
        pool = ThreadPool(3)
        bands = ['g','r','i']
        out = pool.map(self.img, bands)
        pool.close()
        pool.join()
        # Read the trimmed catalog for chip S4 and run on each object to generate mask
        with open("./work/trimcatalog_9999999_4S.pars","r") as f:
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
            print(lines)
            return
            pool = ThreadPool(120)
            out = pool.map(self.mask, lines)
            pool.close()
            pool.join()
        # Remove everything in the work directory
        for f in os.listdir("./work"):
            try: os.remove(os.path.abspath(f))
            except: pass

if __name__ == "__main__":
    for i in range(2,100):
        s = PhoSimSet(i)
        s.simulate()
