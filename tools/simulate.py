import os
import time
import subprocess
import multiprocessing as mp
import numpy as np
from multiprocessing.dummy import Pool as ThreadPool
from astropy.io import fits
from astropy.io.fits import getdata

PHOSIM_DIR = os.path.abspath("../../phosim_release")
TRAIN_DIR = os.path.abspath("../trainingset")
COMMANDFILE_IMG = os.path.abspath("training")
COMMANDFILE_MASK = os.path.abspath("training_nobg")
os.chdir(PHOSIM_DIR)

def bash(command,print_out=True):
    if print_out: print(command)
    try: return subprocess.call(command.split())
    except: return -1

class PhoSimSet:

    def __init__(self,set_num,instrument,bands,exptimes,maglim,train_dir="../trainingset/"):
        self.instrument = instrument
        self.bands = bands
        self.maglim = maglim
        self.set = set_num
        self.seed = 0 + set_num # Use a random scene each time
        self.train_dir = train_dir
        # something roughly like DES foorprint as a test (avoiding galacitic plane)
        self.ra = np.random.uniform(0,60) # deg
        self.dec = np.random.uniform(-70,10) # deg
        # Map each band to a filter number, obsID, and exposure time (WARNING: assuming phosim ugrizy=012345!)
        bands_list = ['u','g','r','i','z','Y']
        self.filterid = dict(zip(bands_list, [0,1,2,3,4,5]))
        self.obsid = dict(zip(bands_list, [9999999,9999998,9999997,9999996,9999995,9999994]))
        self.exptime = dict(zip(self.bands, exptimes))

    def mask(self,line):
        out_dir = os.path.abspath(os.path.join(self.train_dir,"set_%d/" % self.set))
        if len(line.strip()) != 0:
            i = line.split()[1].split('.')[0]
            if line.split()[-3] == "star":
                obj_class = "star"
            else:
                obj_class = "gal"
            # Make a new catalog for each source in the image
            band = self.bands[0] # Use the first band for each mask image
            exptime = 250 # Seconds -- Don't need to use full exposure time
            with open("./examples/obj%s" % i,"w+") as f:
                f.write("rightascension %f\ndeclination %f\nfilter %d\nvistime %f\nnsnap 1\nobshistid %s\nseed %d\n" % (self.ra,self.dec,self.filterid[band],self.exptime[band],i,self.seed))
                f.write(line)
                f.write("\n")
            # Now run PhoSim with this single object and no background to make mask
            bash("./phosim examples/obj%s -c %s -i %s -t %d -e 0" % (i,COMMANDFILE_MASK,self.instrument,int(np.sqrt(mp.cpu_count()))))
            out_to = os.path.abspath(os.path.join(out_dir,"%s%s.fits.gz" % (obj_class,i)))
            out_from = "./output/%s_e_%s_f2_4S_E000.fits.gz" % (self.instrument,i)
            try:
                bash("mv %s %s" % (out_from, out_to))
            except: # off chip
                pass
        print('thread %s done' % i)
        return

    def img(self,band):
        m0 = 15.0
        out_dir = os.path.abspath(os.path.join(self.train_dir,"set_%d/" % self.set))
        out_to = os.path.abspath(os.path.join(out_dir,"img_%s.fits.gz" % band))
        # setup catalog file
        with open("examples/maskrcnn_catalog_%s" % band,"w+") as f:
            print((self.ra,self.dec,self.filterid[band],self.exptime[band],self.obsid[band],self.seed,m0,self.maglim,m0,self.maglim))
            f.write("rightascension %f\ndeclination %f\nfilter %d\nvistime %f\nnsnap 1\nobshistid %d\nseed %d\nstars %f %f 0.1\ngalaxies %f %f 0.1" % (self.ra,self.dec,self.filterid[band],self.exptime[band],self.obsid[band],self.seed,m0,self.maglim,m0,self.maglim))
        # Run PhoSim
        bash("./phosim examples/maskrcnn_catalog_%s -c %s -i %s -t %d -e 0" % (band,COMMANDFILE_IMG,self.instrument,mp.cpu_count()//3))
        # Find chips in output and make set directories for each
        bash("mv ./output/%s_e_%d_f%d_4S_E000.fits.gz %s" % (self.instrument,self.obsid[band],self.filterid[band],out_to))
        return

    def simulate(self):
        out_dir = os.path.abspath(os.path.join(self.train_dir,"set_%d/" % self.set))
        if os.path.exists(out_dir):
            print("Error: Set %d already exists." % self.set)
            return
        os.mkdir(out_dir)
        # Run PhoSim on chip S4 with typical stars and galaxies
        pool = ThreadPool(len(self.bands))
        out = pool.map(self.img, self.bands)
        pool.close()
        pool.join()
        # Read the trimmed catalog for chip S4 and run on each object to generate mask
        with open("./work/trimcatalog_%d_4S.pars" % self.obsid[self.bands[0]],"r") as f:
            lines = []
            readlines = f.readlines()
            for i,line in enumerate(readlines):
                if len(line.strip()) != 0:
                    obj_id = line.split()[1]
                    if ".1" in obj_id: # 1st component of galaxy
                        # add second component
                        lines.append(line+"\n"+readlines[i+1])
                    elif ".2" in obj_id:
                        continue # already added 2nd component
                    else: # star
                        lines.append(line)
            pool = ThreadPool(int(np.sqrt(mp.cpu_count())))
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
        hdul.writeto(os.path.join(out_dir,setdir,"masks.fits"),overwrite=True)
        del hdul

if __name__ == "__main__":
    import argparse
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Simulate training set data using PhoSim.')
    parser.add_argument("instrument",type=str,help="Which PhoSim telescope/instrument to use.")
    parser.add_argument("--bands",type=str,default='grz',help="Which bands to use, e.g. 'grz' or 'grizY'.")
    parser.add_argument("--exptimes",type=str,default='120,120,120',help="Exposure times to use in each band, e.g. '120,120,120'.")
    parser.add_argument("--nset",type=int,default=1,help="Number of sets of images to simulate.")
    parser.add_argument("--maglim",type=float,default='23.0',help="Limiting magnitude.")
    args = parser.parse_args()
    exptimes = np.array(args.exptimes.split(","),dtype=float)

    print('------------------------------------------------------------------------------------------')
    print('Simulate PhoSim Training Data')
    print('------------------------------------------------------------------------------------------')
    print('Instrument: %s' % args.instrument)
    print('Bands: %s' % args.bands)
    print('Exp. Times: %s' % args.exptimes)
    print('Number of sets: %d' % args.nset)
    print('Limiting mag: %2.1f' % args.maglim)
    print('------------------------------------------------------------------------------------------')

    # set loop
    for i in range(args.nset):
        s = PhoSimSet(i,args.instrument,args.bands,exptimes,args.maglim,TRAIN_DIR)
        s.simulate()
    combine_masks(TRAIN_DIR)
