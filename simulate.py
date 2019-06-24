import os
import time
import subprocess
import multiprocessing as mp
import numpy as np
from multiprocessing.dummy import Pool as ThreadPool

PHOSIM_DIR = os.path.abspath("../phosim_core")
os.chdir(PHOSIM_DIR)
TRAIN_DIR = "../deblend_maskrcnn/trainingset"

def bash(command,print_out=True):
    if print_out: print(command)
    try: return subprocess.call(command.split())
    except: return -1

class PhoSimSet:

    def __init__(self,set_num,train_dir="../trainingset/"):
        self.set = set_num
        self.seed = 1000 + set_num
        self.train_dir = train_dir
        # something roughly like DES foorprint as a test
        self.ra = np.random.uniform(0,60) # deg
        self.dec = np.random.uniform(-70,10) # deg

    def mask(self,line):
        out_dir = os.path.abspath(os.path.join(self.train_dir,"set_%d/" % self.set))
        if len(line.strip()) != 0:
            i = line.split()[1].split('.')[0]
            if line.split()[-3] == "star":
                obj_class = "star"
            else:
                obj_class = "gal"
            # Make a new catalog for each source in the image
            with open("./examples/obj%s" % i,"w+") as fi:
                fi.write("rightascension %f\ndeclination %f\nfilter 2\nvistime 150.0\nnsnap 1\nobshistid %s\nseed %d\n" % (self.ra,self.dec,i,self.seed))
                fi.write(line)
                fi.write("\n")
            # Now run PhoSim with this single object and no background to make mask
            bash("./phosim examples/obj%s -c examples/training_nobg -i decam -t 12 -e 0" % i)
            out_to = os.path.abspath(os.path.join(out_dir,"%s%s.fits.gz" % (obj_class,i)))
            out_from = "./output/decam_e_%s_f2_4S_E000.fits.gz" % i
            try:
                bash("mv %s %s" % (out_from, out_to))
            except: # off chip
                pass
        print('thread %s done' % i)
        return

    def img(self,band):
        m0 = 12.0
        m1 = 23.0
        out_dir = os.path.abspath(os.path.join(self.train_dir,"set_%d/" % self.set))
        out_to = os.path.abspath(os.path.join(out_dir,"img_%s.fits.gz" % band))
        with open("examples/maskrcnn_catalog_%s" % band,"w+") as f:
            if band == "g":
                f.write("rightascension %f\ndeclination %f\nfilter 1\nvistime 150.0\nnsnap 1\nobshistid 9999999\nseed %d\nstars %f %f 0.1\ngalaxies %f %f 0.1" % (self.ra,self.dec,self.seed,m0,m1,m0,m1))
            elif band == "r":
                f.write("rightascension %f\ndeclination %f\nfilter 2\nvistime 150.0\nnsnap 1\nobshistid 9999998\nseed %d\nstars %f %f 0.1\ngalaxies %f %f 0.1" % (self.ra,self.dec,self.seed,m0,m1,m0,m1))
            elif band == "z":
                f.write("rightascension %f\ndeclination %f\nfilter 4\nvistime 150.0\nnsnap 1\nobshistid 9999997\nseed %d\nstars %f %f 0.1\ngalaxies %f %f 0.1" % (self.ra,self.dec,self.seed,m0,m1,m0,m1))
        # Run PhoSim
        bash("./phosim examples/maskrcnn_catalog_%s -c examples/training -i decam -t %d -e 0" % (band,mp.cpu_count()//3))
        if band == 'g':
            bash("mv ./output/decam_e_9999999_f1_4S_E000.fits.gz %s" % out_to)
        elif band == 'r':
            bash("mv ./output/decam_e_9999998_f2_4S_E000.fits.gz %s" % out_to)
        elif band == 'z':
            bash("mv ./output/decam_e_9999997_f4_4S_E000.fits.gz %s" % out_to)
        return

    def simulate(self):
        num_cpu = mp.cpu_count()
        if num_cpu > 50:
            num_cpu = 50  # avoid NFS overload/OOM error
        out_dir = os.path.abspath(os.path.join(self.train_dir,"set_%d/" % self.set))
        if os.path.exists(out_dir):
            print("Error: Set %d already exists." % self.set)
            return
        os.mkdir(out_dir)
        # Run PhoSim on chip S4 with typical stars and galaxies
        pool = ThreadPool(3)
        bands = ['g','r','z']
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
                    elif ".2" in obj_id:
                        continue # already added 2nd component
                    else: # star
                        lines.append(line)
            pool = ThreadPool(12)
            out = pool.map(self.mask, lines)
            pool.close()
            pool.join()
        # Remove everything in the work directory
        os.chdir(PHOSIM_DIR)
        bash("rm -r ./examples/obj*")
        bash("rm -r ./work/*.pars")

if __name__ == "__main__":
    # set loop
    for i in range(432,500):
        s = PhoSimSet(i,TRAIN_DIR)
        s.simulate()
