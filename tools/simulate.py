import os
import time
import subprocess
import multiprocessing as mp
from glob import glob as glob
from itertools import repeat
import random
import numpy as np
from multiprocessing.dummy import Pool as ThreadPool
from astropy.io import fits
from astropy.io.fits import getdata

PHOSIM_DIR = os.path.abspath("../../phosim")
TRAIN_DIR = os.path.abspath("./trainingset")
COMMANDFILE_IMG = os.path.abspath("../../phosim/examples/training")
COMMANDFILE_MASK = os.path.abspath("../../phosim/examples/training_nobg")
os.chdir(PHOSIM_DIR)


def bash(command,print_out=True):
    if print_out: print(command)
    try: return subprocess.call(command.split())
    except: return -1

class PhoSimSet:

    def __init__(self,npointing,nset,instrument,fov,bands,exptimes,maglim,train_dir,nproc,seed):
        self.npointing = npointing # Pointing iterator
        self.nset = nset # Set iterator (total pointings + chips)
        self.instrument = instrument
        self.bands = bands
        self.maglim = maglim
        self.seed = seed + npointing # Iterate the random seed for each pointing
        self.train_dir = train_dir
        self.nproc = nproc
        # something roughly like DES foorprint as a test (avoiding galacitic plane)
        self.ra = np.random.uniform(0,60) # deg
        self.dec = np.random.uniform(-70,10) # deg
        self.fov = fov # Degrees to simulate field of stars and galaxies 
        # Map each band to a filter number, obsID, and exposure time
        if self.instrument == 'subaru_hsc': # grizY
            bands_list = ['g','r','i','z','Y']
        else: # Assume ugriz (LSST, decam, etc.)
            bands_list = ['u','g','r','i','z','Y']
        self.filterid = dict(zip(bands_list, range(len(bands_list))))
        self.obsid = dict(zip(bands_list, range(9999999,999999-len(bands_list),-1)))
        self.exptime = dict(zip(self.bands, exptimes))

    def mask(self,line,chip,setnum):
        out_dir = os.path.abspath(os.path.join(self.train_dir,"set_%d/" % setnum))
        if len(line.strip()) != 0:
            # Object number from trimcatalog line gives us the obsID (and therefore the filename)
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
            bash("./phosim examples/obj%s -c %s -i %s -t %d -e 0" % (i,COMMANDFILE_MASK,self.instrument,int(np.sqrt(self.nproc))))
            out_to = os.path.abspath(os.path.join(out_dir,"%s%s.fits.gz" % (obj_class,i)))
            out_from = "./output/%s_e_%s_f%d_%s_E000.fits.gz" % (self.instrument,i,self.filterid[band],chip)
            try:
                bash("mv %s %s" % (out_from, out_to))
            except: # off chip
                pass
        print('thread %s done' % i)
        return

    def img(self,band):
        m0 = 15.0
        # setup catalog file
        with open("examples/maskrcnn_catalog_%s" % band,"w+") as f:
            f.write("rightascension %f\ndeclination %f\nfilter %d\nvistime %f\nnsnap 1\nobshistid %d\nseed %d\nstars %f %f %f\ngalaxies %f %f %f" % (self.ra,self.dec,self.filterid[band],self.exptime[band],self.obsid[band],self.seed,m0,self.maglim,self.fov,m0,self.maglim,self.fov))
        # Run PhoSim
        bash("./phosim examples/maskrcnn_catalog_%s -c %s -i %s -t %d -e 0" % (band,COMMANDFILE_IMG,self.instrument,self.nproc//len(self.bands)))
        # Find chips in output and make set directories for each
        fs = glob("./output/%s_e_%d_f%d_*_E000.fits.gz" % (self.instrument,self.obsid[band],self.filterid[band]))
        chips = [f.split('_')[-3] + '_' + f.split('_')[-2] for f in fs] 
        sets = range(self.nset,self.nset+len(chips))
        self.nset = sets[-1] # Iterate to max
        # Map chip name to set dir
        self.sets = dict(zip(chips,sets))
        # Set (chip) loop
        for i,f in enumerate(fs):
            print('SETS: ', sets)
            # Move and rename final images in each band
            out_dir = os.path.abspath(os.path.join(self.train_dir,"set_%d/" % (sets[i])))
            out_to = os.path.abspath(os.path.join(out_dir,"img_%s.fits.gz" % band))
            os.makedirs(out_dir,exist_ok=True)
            bash("mv %s %s" % (f,out_to))
        return

    def simulate(self):
        # Clean work dir
        print('Cleaning work directory.')
        [os.remove(f) for f in glob(os.path.abspath("./work/*"))]
        print('Simulating.')
        # Run PhoSim pointing with typical stars and galaxies
        pool = ThreadPool(len(self.bands))
        out = pool.map(self.img, self.bands)
        pool.close()
        pool.join()
        # Read the trimmed catalog for chip * and run on each object to generate mask
        fs = glob("./work/trimcatalog_%d_*.pars" % self.obsid[self.bands[0]])
        # Set (chip) loop
        for f in fs:
            chip = f.split('_')[-2] + '_' + f.split('_')[-1].split('.pars')[0]
            # If chip not in sets dict (likely no sources on the chip)
            if not chip in self.sets.keys():
                continue
            setnum = self.sets[chip]
            with open(os.path.abspath(f),"r") as f:
                lines = []
                readlines = f.readlines()
                # Source loop
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
                # Simulate masks
                pool = ThreadPool(int(np.sqrt(self.nproc)))
                out = pool.starmap(self.mask,zip(lines,repeat(chip),repeat(setnum)))
                pool.close()
                pool.join()
        return self.nset

def combine_masks(out_dir):
    # combine masks into one multi-extension HDU
    fs = glob(os.path.abspath(os.path.join(out_dir,'set_*')))
    # Set loop
    for i,f in enumerate(fs):
        hdul = fits.HDUList()
        # mask image loop
        for image in os.listdir(f):
            if (image.endswith('.fits.gz') or image.endswith('.fits')) and not 'img' in image:
                image = os.path.join(f,image)
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
                os.remove(image)
        hdul.writeto(os.path.join(f,"masks.fits"),overwrite=True)
        del hdul

if __name__ == "__main__":
    import argparse
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Simulate training set data using PhoSim.')
    parser.add_argument("instrument",type=str,help="Which PhoSim telescope/instrument to use.")
    parser.add_argument("--fov",type=float,default=0.1,help="Field of view in degrees of catalog to generate, e.g. 0.01 to test.")
    parser.add_argument("--bands",type=str,default='grz',help="Which bands to use, e.g. 'grz' or 'grizY'.")
    parser.add_argument("--exptimes",type=str,default='120,120,120',help="Exposure times in seconds to use in each band, e.g. '120,120,120'.")
    parser.add_argument("--npoint",type=int,default=1,help="Number of sets of images to simulate.")
    parser.add_argument("--maglim",type=float,default='23.0',help="Limiting magnitude.")
    parser.add_argument("--nproc",type=int,default=mp.cpu_count(),help="Number of processes.")
    parser.add_argument("--seed",type=int,default=random.randrange(0,9999999),help="Random number generator seed (default=random).")
    args = parser.parse_args()
    exptimes = np.array(args.exptimes.split(","),dtype=float)
    
    print('------------------------------------------------------------------------------------------')
    print('Simulate PhoSim Training Data')
    print('------------------------------------------------------------------------------------------')
    print('Output directory: %s' % TRAIN_DIR)
    print('Instrument: %s' % args.instrument)
    print('FoV: %1.2f' % args.fov)
    print('Bands: %s' % args.bands)
    print('Exp. Times: %s' % args.exptimes)
    print('Number of pointings: %d' % args.npoint)
    print('Limiting mag: %2.1f' % args.maglim)
    print('Number of processes: %d' % args.nproc)
    print('Seed: %d' % args.seed)
    print('------------------------------------------------------------------------------------------')

    # WARNING: Be sure to remove lines 740-742 in phosim's phosim.py to keep the trimcatalog file from being removed!

    random.seed(args.seed)

    nset = 0 # Total set number counter
    # Pointings loop
    for i in range(args.npoint):
        s = PhoSimSet(i,nset,args.instrument,args.fov,args.bands,exptimes,args.maglim,TRAIN_DIR,args.nproc,args.seed)
        nset = s.simulate()
    combine_masks(TRAIN_DIR)
