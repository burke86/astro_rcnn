import os
import time
import subprocess
import multiprocessing as mp
import numpy as np
from multiprocessing.dummy import Pool as ThreadPool
from astropy.io import fits
from astropy.io.fits import getdata

PHOSIM_DIR = os.path.abspath("/Users/anshul/Downloads/phosim")
TRAIN_DIR = os.path.abspath("trainingset")
os.chdir(PHOSIM_DIR)

def bash(command,print_out=True):
    if print_out: print(command)
    try: return subprocess.call(command.split())
    except: return -1

class PhoSimSet:

    def __init__(self,set_num,train_dir="../trainingset/"):
        self.set = set_num
        self.seed = 0 + set_num
        self.train_dir = train_dir
        # something roughly like DES foorprint as a test (Have to change for hsc?)
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
                fi.write("rightascension %f\ndeclination %f\nfilter 2\nvistime 120.\nnsnap 1\nobshistid %s\nseed %d\n" % (self.ra,self.dec,i,self.seed))
                fi.write(line)
                fi.write("\n")
            # Now run PhoSim with this single object and no background to make mask
            bash("./phosim examples/obj%s -c examples/training_nobg -i subaru_hsc -s 1_06 -t 12 -e 0" % i)
            out_to = os.path.abspath(os.path.join(out_dir,"%s%s.fits.gz" % (obj_class,i)))
            out_from = "./output/decam_e_%s_f2_4S_E000.fits.gz" % i
            try:
                bash("mv %s %s" % (out_from, out_to))
            except: # off chip
                pass
        print('thread %s done' % i)
        return

    def img(self,band):
        m0 = 15.0
        m1 = 28.0
        out_dir = os.path.abspath(os.path.join(self.train_dir,"set_%d/" % self.set))
        out_to = os.path.abspath(os.path.join(out_dir,"img_%s.fits.gz" % band))
        with open("examples/maskrcnn_catalog_%s" % band,"w+") as f:
            if band == "g":
                f.write("rightascension %f\ndeclination %f\nfilter 1\nvistime 120.\nnsnap 1\nobshistid 9999999\nseed %d\nstars %f %f 0.1\ngalaxies %f %f 0.1" % (self.ra,self.dec,self.seed,m0,m1,m0,m1))
            elif band == "r":
                f.write("rightascension %f\ndeclination %f\nfilter 2\nvistime 120.\nnsnap 1\nobshistid 9999998\nseed %d\nstars %f %f 0.1\ngalaxies %f %f 0.1" % (self.ra,self.dec,self.seed,m0,m1,m0,m1))
            elif band == "z":
                f.write("rightascension %f\ndeclination %f\nfilter 4\nvistime 120.\nnsnap 1\nobshistid 9999997\nseed %d\nstars %f %f 0.1\ngalaxies %f %f 0.1" % (self.ra,self.dec,self.seed,m0,m1,m0,m1))
        # Run PhoSim
        bash("./phosim examples/maskrcnn_catalog_%s -c examples/training -i subaru_hsc -s 1_06 -t %d -e 0" % (band,mp.cpu_count()//3))
        print("DONE")
        return
        if band == 'g':
            bash("mv ./output/lsst_a_0000_f4_R22_S11_C02_E000.fits.gz %s" % out_to)
        elif band == 'r':
            bash("mv ./output/lsst_a_0000_f4_R22_S11_C01_E000.fits.gz %s" % out_to)
        elif band == 'z':
            bash("mv ./output/lsst_a_0000_f4_R22_S11_C00_E000.fits.gz %s" % out_to)
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
        return
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

def combine_masks(out_dir):
    # combine masks into one multi-extension HDU
    num_sets = 0
    for setdir in os.listdir(out_dir):
        if 'set_' in setdir:
            num_sets += 1

    for i in range(num_sets):
        hdul = fits.HDUList()
        setdir = 'set_%d' % i
        for image in os.listdir(os.path.join(out_dir,setdir)):
            if (image.endswith('.fits.gz') or image.endswith('.fits')) and not 'img' in image:
                image = os.path.join(out_dir,setdir,image)
                data = getdata(image)
                # all zeros
                if not np.any(data): continue
                hdr = fits.Header()
                if 'star' in image:
                    hdr["CLASS_ID"] = 1
                elif "gal" in image:
                    hdr["CLASS_ID"] = 2
                else: continue
                hdr["BITPIX"] = 8
                hdul.append(fits.ImageHDU(data,header=hdr))
                del data
        print(os.path.join(out_dir,setdir,"masks.fits"))
        hdul.writeto(os.path.join(out_dir,setdir,"masks.fits"),overwrite=True)
        del hdul

if __name__ == "__main__":
    # set loop
    for i in range(0,1):
        s = PhoSimSet(i,TRAIN_DIR)
        s.simulate()
    combine_masks(TRAIN_DIR)
