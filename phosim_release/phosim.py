#!/usr/bin/env python
##
## @package phosim.py
## @file phosim.py
## @brief python methods called by phosim script
##
## @brief Created by:
## @author John R. Peterson (Purdue)
##
## @brief Modified by:
## @author Emily Grace (Purdue)
## @author Nathan Todd (Purdue)
## @author En-Hsin Peng (Purdue)
## @author Glenn Sembroski (Purdue)
## @author Jim Chiang (SLAC)
## @author Jeff Gardner (Google)
## @author Colin Burke (Purdue)
##
## @warning This code is not fully validated
## and not ready for full release.  Please
## treat results with caution.
##
## The intention is that script is the only file necessary to run phosim for
## any purpose-- ranging from a single star on your laptop to full fields
## on large-scale computing clusters.  There is no physics in this
## script.  Its only purpose is to move files to the correct place.
## This file is called by the phosim script

import os
import subprocess
import sys, glob, optparse, shutil
import multiprocessing
import math

## print the usage
def usage():
     script=os.path.abspath('phosim')
     os.system(script+' x --help')

##jobChip is a function that run an individual chip for a single exposure
def jobChip(observationID, cid, eid, filt, nframes, nskip, ngroups, outputDir, binDir, instrDir, instrument='lsst', run_e2adc=True, run_ds9=False):
    fid = observationID + '_' + cid + '_' + eid
    segfile = instrDir+'/segmentation.txt'
    runProgram("raytrace < raytrace_"+fid+".pars", binDir)
    removeFile('raytrace_'+fid+'.pars')
    nsamples = nframes + nskip
    # try to find first image before we run e2adc program
    if (nsamples > 1 or ngroups > 1): # sequence mode
        nstop = ngroups + 1 # include reference frame as "group 0"
        eImage = instrument+'_e_'+observationID+'_f'+filt+'_'+cid+'_'+eid+'_R000.fits.gz'
    else: # frame mode
        nstop = 1
        eImage = instrument+'_e_'+observationID+'_f'+filt+'_'+cid+'_'+eid+'.fits.gz'
    if os.path.exists(eImage):
         if run_e2adc:
              runProgram("e2adc < e2adc_"+fid+".pars", binDir)
              removeFile('e2adc_'+fid+'.pars')
              # move amplifier images to output directory
              for line in open(segfile):
                   aid = line.split()[0]
                   if cid in line and aid != cid:
                        for ngroup in range(nstop):
                             if (nstop == 1):
                                  rawImage = instrument+'_a_'+observationID+'_f'+filt+'_'+aid+'_'+eid+'.fits.gz'
                             else:
                                  rawImage = instrument+'_a_'+observationID+'_f'+filt+'_'+aid+'_'+eid+'_G'+str(ngroup).zfill(3)+'.fits.gz'
                             if os.path.exists(rawImage):
                                  shutil.move(rawImage,outputDir+'/'+rawImage)
         # move electron images to output directory
         nframei = 0
         for nframe in range(ngroups*nsamples - nskip):
              if (nstop != 1):
                  eImage=instrument+'_e_'+observationID+'_f'+filt+'_'+cid+'_'+eid+'_R'+str(nframe).zfill(3)+'.fits.gz'
              if os.path.exists(eImage):
                  shutil.move(eImage, outputDir+'/'+eImage)
              # do not simulate frame skips
              if nframei % nframes == 0:
                  nframe += nskip
                  framei = 0
         if run_ds9:
             fullpath = '/Applications:' + os.environ["PATH"]
             for path in fullpath.split(os.pathsep):
                  if os.path.exists(os.path.join(path, "ds9")):
                       commandName = os.path.join(path, 'ds9 -scale log ' + outputDir + '/' + eImage + ' &')
                       if subprocess.call(commandName, shell=True) != 0:
                            sys.exit(1)

## runProgram function calls each of the phosim programs using subprocess.call
#  it raises an exception and aborts if the return code is non-zero.
def runProgram(command, binDir=None, argstring=None):
    myCommand = command
    if binDir is not None:
        myCommand = os.path.join(binDir, command)
    if argstring is not None:
        myCommand += argstring
    if subprocess.call(myCommand, shell=True) != 0:
        sys.exit(1)

## removeFile deletes files (if they do not exist, it will catch the OSError exception and silently proceed.)
def removeFile(filename):
    try:
         os.remove(filename)
    except OSError:
         pass

## assignPath figures out the whole path
def assignPath(path,phosimDir):
    if path == 'none':
        return 'none'
    full_path = phosimDir + '/' + path
    if os.path.exists(full_path):
        return full_path
    elif os.path.exists(path):
        return os.path.abspath(path)
    raise RuntimeError('%s does not exist.' % path)

## checkPaths makes sure the paths make sense
def checkPaths(opt,phosimDir):
    for x in ['outputDir','workDir','binDir','dataDir','sedDir','imageDir','extraCommands']:
        exec('opt.%s=assignPath(opt.%s,phosimDir)' % (x,x))

    #Now find the instrument directory
    instrDir = os.path.join(opt.dataDir, opt.instrument)
    if not os.path.exists(instrDir):
        raise RuntimeError('The instrument directory %s does not exist.' % instrDir)
    opt.instrDir=instrDir
    opt.instrument = instrDir.split("/")[-1]
    if len(opt.instrument)==0:
        opt.instrument = instrDir.split("/")[-2]


## PhosimFocalplane is a class for handling phosim files and directories for one focalplane.
#    In order to perform all of the preprocessing steps, use doPreproc().  The
#    raytrace steps can then be scheduled (exactly how depends on the 'grid' option)
#    using scheduleRaytrace().  Once the computation is finished, intermediate files
#    can be removed using cleanup().  Thus, an example workflow would be as follows:
#    focalplane = PhosimFocalplane(phosimDir, outputDir, workDir, binDir, dataDir, instrDir, opt.grid, grid_opts)
#        focalplane.doPreproc(instanceCatalog, extraCommands, sensor)
#        focalplane.scheduleRaytrace(instrument, run_e2adc, keep_screens)
#        focalplane.cleanup(keep_screens)
class PhosimFocalplane(object):

     ## Constructor for PhosimFocalplane
     #  grid:      'no', 'condor', 'cluster'
     #  grid_opts: A dictionary to supply grid options.  Exactly which options
     #  depends on the value of 'grid':
     #  'no':      'numproc' = Number of threads used to execute raytrace.
     #  'condor':  'universe' = Condor universe ('vanilla', 'standard', etc)
     #  'cluster': 'script_writer' = callback to generate raytrace batch scripts
     #             'submitter' = optional callback to submit the job
     def __init__(self, phosimDir, opt, grid_opts={}, visualize=False):

        self.phosimDir = phosimDir
        self.outputDir = opt.outputDir
        self.workDir = opt.workDir
        self.binDir = opt.binDir
        self.dataDir = opt.dataDir
        self.instrDir = opt.instrDir
        self.sedDir = opt.sedDir
        self.imageDir = opt.imageDir
        self.numthread = opt.numthread
        self.flatdir = False
        self.tarfile = False
        self.extraCommands = None
        self.instanceCatalog = None
        self.userCatalog = None
        self.chipID = None
        self.runFlag = None
        self.devmaterial = None
        self.devtype = None
        self.devvalue = None
        self.devmode = None
        self.grid = opt.grid
        self.grid_opts = grid_opts
        self.execEnvironmentInitialized = False
        self.eventfile = visualize
        if self.grid == 'condor':
            assert 'universe' in self.grid_opts
            self.flatdir=True if self.grid_opts['universe'] == 'vanilla' else False
        elif self.grid == 'diagrid':
            self.flatdir=True
            self.tarfile=True

     ## doPreproc is a method to run all of the non-chip steps.
     def doPreproc(self, instanceCatalog, extraCommands, sensor):
        self.loadInstanceCatalog(instanceCatalog, extraCommands)
        os.chdir(self.workDir)
        self.writeInputParamsAndCatalogs()
        self.generateAtmosphere()
        self.generateInstrumentConfig()
        self.trimObjects(sensor)

     ## Parses the instance catalog
     def loadInstanceCatalog(self, instanceCatalog, extraCommands):
        self.instanceCatalog = instanceCatalog
        self.extraCommands = extraCommands
        self.nframes = 1
        self.nskip = 0
        defaultCatalog=open(os.path.join(self.phosimDir,'default_instcat')).readlines()
        self.userCatalog=open(instanceCatalog).readlines()
        moveString = ''
        self.catgen = 0 # flag indicating we are to generate some object catalogs
        for line in defaultCatalog+self.userCatalog:
             lstr=line.split()
             if "obshistid" in line:
                  self.observationID=lstr[1]
             elif "moonra" in line:
                  self.moonra=lstr[1]
             elif "moondec" in line:
                  self.moondec=lstr[1]
             elif "sunalt" in line:
                  self.solaralt=lstr[1]
             elif "moonalt" in line:
                  self.moonalt=lstr[1]
             elif "dist2moon" in line:
                  self.moondist=lstr[1]
             elif "moonphase" in line:
                  self.phaseang=lstr[1]
             elif "mjd" in line:
                  self.tai=lstr[1]
             elif "seeing" in line:
                  self.constrainseeing=lstr[1]
             elif "airglow" in line:
                  self.constrainairglow=lstr[1]
             elif "clouds" in line:
                  self.constrainclouds=lstr[1]
             elif "rottelpos" in line:
                  self.spiderangle=lstr[1]
             elif "Azimuth" in line or "azimuth" in line:
                  self.azimuth=lstr[1]
             elif "Altitude" in line or "altitude" in line:
                  self.altitude=lstr[1]
             elif "rotskypos" in line:
                  self.rotationangle=lstr[1]
             elif "Unrefracted_RA" in line or "rightascension" in line:
                  self.pointingra=lstr[1]
             elif "Unrefracted_Dec" in line or "declination" in line:
                  self.pointingdec=lstr[1]
             elif "SEED" in line  or "seed" in line:
                  self.obsseed=lstr[1]
             elif "date" in line:
                  self.monthnum=lstr[1].split('/')[1]
             elif "filter" in line:
                  self.filt=lstr[1]
             elif "VISTIME" in line or "vistime" in line:
                  self.vistime=float(lstr[1])
             elif "NSNAP" in line or "nsnap" in line:
                  self.nsnap=int(float(lstr[1]))
             elif "NFRAMES" in line or "nframes" in line:
                  self.nframes=int(float(lstr[1]))
             elif "NSKIP" in line or "nskip" in line:
                  self.nskip=int(float(lstr[1]))
             elif "MINSOURCE" in line or "minsource" in line:
                  self.minNumSources=int(float(lstr[1]))
             elif "CAMCONFIG" in line or "camconfig" in line:
                  self.camconfig=int(float(lstr[1]))
             elif "DOMEINT" in line or "domeint" in line:
                  self.domeint=float(lstr[1])
             elif "DOMEWAV" in line or "domewav" in line:
                  self.domewav=float(lstr[1])
             elif "TELCONFIG" in line or "telconfig" in line:
                  self.telconfig=int(float(lstr[1]))
             elif "TEMPERATURE" in line or "temperature" in line:
                  self.temperature=lstr[1]
             elif "TEMPVAR" in line or "tempvar" in line:
                  self.tempvar=lstr[1]
             elif "PRESSURE" in line or "pressure" in line:
                  self.pressure=lstr[1]
             elif "OVERDEPBIAS" in line or "overdepbias" in line:
                  self.overdepbias=lstr[1]
             elif "CONTROL" in line or "control" in line:
                  self.control=lstr[1]
             elif "move" in line:
                  moveString += line.strip().split("move")[1]
             elif "focus" in line:
                  self.focus = line.strip().split("focus")[1]
             elif "stars" in line:
                  self.catgenStars = line.strip().split("stars")[1]
                  self.catgen=1
             elif "stargrid" in line:
                  self.catgenStarGrid =  line.strip().split("stargrid")[1]
                  self.catgen=1
             elif "galaxies" in line:
                  self.catgenGalaxies =  line.strip().split("galaxies")[1]
                  self.catgen=1
        self.actuator=moveString + '\n'
        self.throughputfile=0
        self.centroidfile=0
        self.opdfile=0
        if extraCommands != 'none':
             # look for relevant physics commands
             for line in open(extraCommands):
                  lstr=line.split()
                  if "extraid" in line:
                       self.extraid=lstr[1]
                       self.observationID=self.observationID+self.extraid
                  if "eventfile" in line:
                       self.eventfile=int(float(lstr[1]))
                  if "throughputfile" in line:
                       self.throughputfile=int(float(lstr[1]))
                  if "centroidfile" in line:
                       self.centroidfile=int(float(lstr[1]))
                  if "opd" in line:
                       self.opdfile=int(float(lstr[1]))
             # catch problem with appending pars file without endline character in command file
             if line[-1] != os.linesep:
                  print('Error: No endline character found in command file.')
                  sys.exit()
 
     ## writeInputParamsAndCatalogs encapsulate the two instance catalog processing functions.
     def writeInputParamsAndCatalogs(self):
        self.writeInputParams()
        self.writeCatalogList()

     ## writeInputParams takes some of the parsed input parameters out of the
     # instance catalog and puts them in a single file.  We mainly have to read this
     # in and write out the same file because the instance catalog format cannot change
     #rapidly due to compatibility with opSim and catSim.
     def writeInputParams(self):
        self.inputParams='obs_'+self.observationID+'.pars'
        pfile=open(self.inputParams,'w')
        pfile.write("obshistid %s\n" % self.observationID)
        if hasattr(self, 'moonra'):
             pfile.write("moonra %s\n" % self.moonra)
        if hasattr(self, 'moondec'):
             pfile.write("moondec %s\n" % self.moondec)
        if hasattr(self, 'solaralt'):
             pfile.write("solaralt %s\n" % self.solaralt)
        if hasattr(self, 'moonalt'):
             pfile.write("moonalt %s\n" % self.moonalt)
        if hasattr(self, 'moondist'):
             pfile.write("moondist %s\n" % self.moondist)
        if hasattr(self, 'phaseang'):
             pfile.write("phaseang %s\n" % self.phaseang)
        if hasattr(self, 'tai'):
             pfile.write("tai %s\n" % self.tai)
        if hasattr(self, 'azimuth'):
             pfile.write("azimuth %s\n" % self.azimuth)
        if hasattr(self, 'altitude'):
             pfile.write("altitude %s\n" % self.altitude)
        if hasattr(self, 'pointingra'):
             pfile.write("pointingra %s\n" % self.pointingra)
        if hasattr(self, 'pointingdec'):
             pfile.write("pointingdec %s\n" % self.pointingdec)
        if hasattr(self, 'monthnum'):
             pfile.write("monthnum %s\n" % self.monthnum)
        if hasattr(self, 'spiderangle'):
             pfile.write("spiderangle %s\n" % self.spiderangle)
        if hasattr(self, 'rotationangle'):
             pfile.write("rotationangle %s\n" % self.rotationangle)
        if hasattr(self, 'filt'):
             pfile.write("filter %s\n" % self.filt)
        if hasattr(self, 'catgenStars'):
             pfile.write("stars %s\n" % self.catgenStars)
        if hasattr(self, 'catgenStarGrid'):
             pfile.write("stargrid %s\n" % self.catgenStarGrid)
        if hasattr(self, 'catgenGalaxies'):
             pfile.write("galaxies %s\n" % self.catgenGalaxies)
        if hasattr(self, 'temperature'):
             pfile.write("temperature %s\n" % self.temperature)
        if hasattr(self, 'pressure'):
             pfile.write("pressure %s\n" % self.pressure)
        if hasattr(self, 'tempvar'):
             pfile.write("tempvar %s\n" % self.tempvar)
        if hasattr(self, 'focus'):
             pfile.write("focus %s\n" % self.focus)

        pfile.write("constrainseeing %s\n" % self.constrainseeing)
        pfile.write("constrainairglow %s\n" % self.constrainairglow)
        pfile.write("constrainclouds %s\n" % self.constrainclouds)
        pfile.write("overdepbias %s\n" % self.overdepbias)
        pfile.write("move %s" % self.actuator)
        pfile.write("control %s\n" % self.control)
        pfile.write("obsseed %s\n" % self.obsseed)
        pfile.write("vistime %g\n" % self.vistime)
        pfile.write("camconfig %s\n"% self.camconfig)
        pfile.write("outputdir %s\n" % self.outputDir)
        pfile.write("seddir %s\n" % self.sedDir)
        pfile.write("imagedir %s\n" % self.imageDir)
        pfile.write("datadir %s\n" % self.dataDir)
        pfile.write("instrdir %s\n" % self.instrDir)
        pfile.write("bindir %s\n" % self.binDir)
        pfile.write("thread %d\n" % self.numthread)
        pfile.write("telconfig %d\n" % self.telconfig)
        pfile.write("minsource %d\n" % self.minNumSources)
        pfile.write("domelight %g\n" % self.domeint)
        pfile.write("domewave %g\n" % self.domewav)
        if self.flatdir:
             pfile.write("flatdir 1\n")
        if self.tarfile:
             pfile.write("tarfile 1\n")
        if self.eventfile:
             pfile.write("eventfile 1\n")
        pfile.close()

     ## writeCatalogList simply makes a list of possible
     # sub-catalogs (using the includeobj option) or lists of
     # objects put in the instance catalog.  The former is useful
     # for 1000s of objects, whereas the latter is useful for entire
     # focalplanes (millions).  Hence we support both of these options.
     def writeCatalogList(self):
        assert self.instanceCatalog
        assert self.userCatalog
        l=0
        objectCatalog=open('objectcatalog_'+self.observationID+'.pars','w')
        for line in self.userCatalog:
             if "object" in line:
                  objectCatalog.write(line)
                  l+=1
        objectCatalog.close()
        #OPD
        l2=0
        opdCatalog=open('opdcatalog_'+self.observationID+'.pars','w')
        for line in self.userCatalog:
             if "opd" in line:
                  opdCatalog.write(line)
                  l2+=1
        opdCatalog.close()
        ncat=0
        catalogList=open('catlist_'+self.observationID+'.pars','w')
        if l>0:
             catalogList.write("catalog %d objectcatalog_%s.pars\n" % (ncat,self.observationID))
             ncat=1
        else:
             removeFile('objectcatalog_'+self.observationID+'.pars')
        #opd
        if l2>0:
             catalogList.write("catalogopd opdcatalog_%s.pars\n" % (self.observationID))
             ncat+=1
        else:
             removeFile('opdcatalog_'+self.observationID+'.pars')
        catDir = os.path.dirname(self.instanceCatalog)
        for line in self.userCatalog:
            if "includeobj" in line:
                path = os.path.join(catDir, line.split()[1])
                if not os.path.isabs(catDir):
                    path = os.path.join("..", path)
                catalogList.write("catalog %d %s\n" % (ncat, path))
                ncat+=1
        catalogList.close()

     ## generateAtmosphere runs the atmosphere program
     def generateAtmosphere(self):
          assert self.inputParams
          assert self.extraCommands
          inputParams='obsExtra_'+self.observationID+'.pars'
          pfile=open(inputParams,'w')
          pfile.write(open(self.inputParams).read())
          if self.extraCommands!='none':
               pfile.write(open(self.extraCommands).read())
          pfile.close()
          runProgram("atmosphere < "+inputParams, self.binDir)
          removeFile(inputParams)
          obsFile=('obs_'+self.observationID+'.pars')
          pfile=open(obsFile,'r')
          for line in pfile:
              lstr=line.split()
              if "filter" in line:
                  self.filt=lstr[1]


     ## generateInstrumentConfig runs the instrument program
     def generateInstrumentConfig(self):
          assert self.inputParams
          assert self.extraCommands
          inputParams='obsExtra_'+self.observationID+'.pars'
          pfile=open(inputParams,'w')
          pfile.write(open(self.inputParams).read())
          if self.extraCommands!='none':
               pfile.write(open(self.extraCommands).read())
          pfile.close()
          runProgram("instrument < "+inputParams, self.binDir)
          removeFile(inputParams)

     ## Run a program that will generate user defined simulated catalogs.
     # This can include random stars, and/or  a grid of stars and/or random
     # galaxies. Options for this are specified in the input command file
     def generateSimulatedCatalogs(self):
         # Have we been requested to make simulated catalogs?
         if self.catgen != 0 :
              # All commnad and info need to make catalogs are contained
              #in the obs_*.pars file.  Now run it
              obsFile=('obs_'+self.observationID+'.pars')
              #print("about to run phosimcatgen")
              runProgram("phosimcatgen < " + obsFile, self.binDir)
               
              # We shoud now have a catalog file. Append it to the end of the
              # catlist file. First generate the ncat number of this
              # additional (or perhps only)  catalog
              catgenFile = 'catgen_'+ self.observationID +'.cat'
              
              # We now need to make sure the catlist file exists (it should
              # I think,  but it may be empty
              catListFile='catlist_'+self.observationID+'.pars'
              if os.path.exists(catListFile):
                   if os.path.getsize(catListFile) > 0:
                        catalogList=open(catListFile)
                        lineList = catalogList.readlines()
                        if "catalog" == lineList[-1].split()[0] :
                             ncat = int(lineList[-1].split()[1])
                             ncat+=1
                        catalogList.close()
                   else:
                        ncat=0
              else:
                   ncat=0   #This is just in case

              # Append new catalog to existing catlist file or to a newly
              # created one
              catalogList=open(catListFile,'a+')
              catalogList.write("catalog %d %s\n" % (ncat, catgenFile))
              catalogList.close()
              # And we are done. 

     ## trimObjects runs the trim program
     #  Note this is overly complicated because we want to allow the trimming
     #  on grid computing to be done in groups to reduce the I/O of sending
     #  the entire instance catalog for every chip.  This complex looping
     #  isn't necessary for desktops.
     def trimObjects(self, sensors):

          self.initExecutionEnvironment()

          camstr="%03d" % int(float(bin(self.camconfig).split('b')[1]))
          if self.camconfig==0:
               camstr='111'
          fp=open(self.instrDir+"/focalplanelayout.txt").readlines()
          chipID=[]
          runFlag=[]
          devmaterial=[]
          devtype=[]
          devmode=[]
          devvalue=[]

          #Go through the focalplanelayout.txt filling up the arrays
          for line in fp:
               lstr=line.split()
               addFlag=0
               if "Group0" in line and camstr[2]=='1': addFlag=1
               elif "Group1" in line and camstr[1]=='1': addFlag=1
               elif "Group2" in line and camstr[0]=='1': addFlag=1
               if addFlag==1:
                    chipID.append(lstr[0])
                    runFlag.append(1)
                    devmaterial.append(lstr[6])
                    devtype.append(lstr[7])
                    devmode.append(lstr[8])
                    devvalue.append(float(lstr[9]))

          # OPD
          chipID.append('opd')
          runFlag.append(1)
          devmaterial.append('silicon')
          devtype.append('CCD')
          devmode.append('frame')
          devvalue.append(0)

          # See if we limit ourselves to a specific set of chipID (seperated by "|").
          if sensors != 'all':
               lstr = sensors.split('|')
               for i in range(len(chipID)): runFlag[i]=0
               for j in range(len(lstr)):
                    for i in range(len(chipID)):
                         if lstr[j]==chipID[i]:
                              runFlag[i]=1
                              break

          # OPD
          secondlastchip=chipID[-2]
          lastchip=chipID[-1]
          chipcounter1=0
          chipcounter2=0
          tc=0
          i=0
          trimJobID=[]
          for cid in chipID:
               if chipcounter1==0:
                    jobName='trim_'+self.observationID+'_'+str(tc)
                    inputParams=jobName+'.pars'
                    pfile=open(inputParams,'w')

               pfile.write('chipid %d %s\n' % (chipcounter1,cid))
               chipcounter1+=1
               if runFlag[i]==1:
                    chipcounter2+=1
          #OPD
               if chipcounter1==9 or cid==lastchip or cid==secondlastchip:   #Do groups of 9 to reduce grid computing I/O
                    trimJobID.append('none')
                    pfile.write(open('obs_'+self.observationID+'.pars').read())
                    if self.flatdir:
                         for line in open('catlist_'+self.observationID+'.pars'):
                              lstr=line.split()
                              pfile.write('%s %s %s\n' % (lstr[0],lstr[1],lstr[2].split('/')[-1]))
                    else:
                         pfile.write(open('catlist_'+self.observationID+'.pars').read())
                    pfile.close()
                    if chipcounter2>0:
                         if self.grid in ['no', 'cluster']:
                              runProgram("trim < "+inputParams, self.binDir)
                         elif self.grid == 'condor':
                              nexp=self.nsnap if devtype[i]=='CCD' else int(self.vistime/devvalue[i])
                              condor.writeTrimDag(self,jobName,tc,nexp)
                         elif self.grid == 'diagrid':
                              nexp=self.nsnap if devtype[i]=='CCD' else int(self.vistime/devvalue[i])
                              trimJobID[tc]=diagrid.writeTrimDag(self,jobName,tc,nexp)
                         else:
                              sys.stderr.write('Unknown grid type: %s' % self.grid)
                              sys.exit(-1)
                    if self.grid in ['no', 'cluster'] or (self.grid in ['condor','diagrid'] and chipcounter2==0):
                         removeFile(inputParams)
                    chipcounter1=0
                    chipcounter2=0
                    tc+=1
               i=i+1
          self.chipID = chipID
          self.runFlag = runFlag
          self.devmaterial = devmaterial
          self.devtype = devtype
          self.devvalue = devvalue
          self.devmode = devmode
          self.trimJobID = trimJobID

     #scheduleRaytrace sets up the raytrace & e2adc jobs and also figures out the
     #numbers of exposures to perform.
     def scheduleRaytrace(self, instrument='lsst', run_e2adc=True, keep_screens=False, run_ds9=False, run_phosimvisualizer=False):
        assert self.extraCommands
        chipcounter1=0
        tc=0
        counter=0
        jobs=[]
        rImageArg='' # list of eventfile paths for phosim visualizer
        i=0
        seg=open(self.instrDir+'/segmentation.txt').readlines()
        observationID = self.observationID

        for cid in self.chipID:
            if self.runFlag[i]==1:
                numSources=self.minNumSources
                if self.grid in ['no', 'cluster']:
                    numSources=len(open('trimcatalog_'+observationID+'_'+cid+'.pars').readlines())
                    numSources=numSources-2
                if numSources>=self.minNumSources:
                    nexp=self.nsnap
                    if self.devtype[i]=='CMOS' and self.devmode[i]=='frame':
                        nexp=int(self.vistime/self.devvalue[i])
                    ex=0
                    while ex<nexp:
                        eid="E%03d" % (ex)
                        fid=observationID + '_' + cid + '_' + eid
                        pfile=open('image_'+fid+'.pars','w')
                        pfile.write("chipid %s\n" % cid)
                        pfile.write("exposureid %d\n" % ex)
                        pfile.write("nsnap %d\n" % nexp)
                        pfile.write("nframes %d\n" % self.nframes)
                        pfile.write("nskip %d\n" % self.nskip)
                        if self.nframes == 1 and self.nskip == 0: # override
                          self.devmode[i] = 'frame'
                        if self.devmode[i] == 'frame':
                          self.ngroups = 1
                          self.nframes = 1
                          self.skip = 0
                          nsamples = self.nframes + self.nskip
                        else:
                          nsamples = self.nframes + self.nskip
                          self.ngroups = int(math.ceil((self.vistime + self.devvalue[i]*self.nskip)/(self.devvalue[i]*nsamples)))
                        pfile.close()
                        
                        # PHOTON RAYTRACE
                        pfile=open('raytrace_'+fid+'.pars','w')
                        pfile.write(open('obs_'+observationID+'.pars').read())
                        if os.path.exists('atmosphere_'+observationID+'.pars'):
                             pfile.write(open('atmosphere_'+observationID+'.pars').read())
                        pfile.write(open('optics_'+observationID+'.pars').read())
                        if (cid != 'opd'):
                             pfile.write(open('chip_'+observationID+'_'+cid+'.pars').read())
                        pfile.write(open('image_'+fid+'.pars').read())
                        if self.extraCommands!='none':
                            pfile.write(open(self.extraCommands).read())
                        if self.grid in ['no', 'cluster']:
                            pfile.write(open('trimcatalog_'+observationID+'_'+cid+'.pars').read())
                        pfile.close()

                        # ELECTRONS TO ADC CONVERTER
                        if run_e2adc:
                            pfile=open('e2adc_'+fid+'.pars','w')
                            pfile.write(open('obs_'+observationID+'.pars').read())
                            if (cid != 'opd'):
                                 pfile.write(open('readout_'+observationID+'_'+cid+'.pars').read())
                            if self.extraCommands!='none':
                                 pfile.write(open(self.extraCommands).read())
                            pfile.write(open('image_'+fid+'.pars').read())
                            pfile.close()

                        if self.grid == 'no':
                            p=multiprocessing.Process(target=jobChip,
                                                      args=(observationID,cid,eid,self.filt, self.nframes, self.nskip, self.ngroups, self.outputDir,
                                                            self.binDir, self.instrDir),
                                                      kwargs={'instrument': instrument, 'run_e2adc': run_e2adc, 'run_ds9': run_ds9})
                            jobs.append(p)
                            p.start()
                            counter+=1
                            if counter==self.grid_opts.get('numproc', 1):
                                for p in jobs:
                                    p.join()
                                counter=0
                                jobs=[]
                        elif self.grid == 'cluster':
                            if self.grid_opts.get('script_writer', None):
                                self.grid_opts['script_writer'](observationID, cid, eid, self.filt,
                                                                   self.outputDir, self.binDir, self.dataDir)
                            else:
                                sys.stderr.write('WARNING: No script_writer callback in grid_opts for grid "cluster".\n')
                            if self.grid_opts.get('submitter', None):
                                self.grid_opts['submitter'](observationID, cid, eid)
                            else:
                                sys.stdout.write('No submitter callback in self.grid_opts for grid "cluster".\n')
                        elif self.grid == 'condor':
                            condor.writeRaytraceDag(self,cid,eid,tc,run_e2adc)
                        elif self.grid == 'diagrid':
                            diagrid.writeRaytraceDag(self,cid,eid,tc,run_e2adc)

                        if run_phosimvisualizer:
                            rImage = os.path.join(self.outputDir, instrument+'_r_'+observationID+'_f'+self.filt+'_'+cid+'_'+eid+'.fits')
                            rImageArg += rImage+' '

                        removeFile('image_'+fid+'.pars')
                        ex+=1

            chipcounter1+=1
            if chipcounter1==9:
                tc+=1
                chipcounter1=0

            if self.grid in ['no', 'cluster']:
                if os.path.exists('trimcatalog_'+observationID+'_'+cid+'.pars'):
                    removeFile('trimcatalog_'+observationID+'_'+cid+'.pars')
            removeFile('readout_'+observationID+'_'+cid+'.pars')
            removeFile('chip_'+observationID+'_'+cid+'.pars')
            i+=1

        removeFile('obs_'+observationID+'.pars')
        if not keep_screens:
             removeFile('atmosphere_'+observationID+'.pars')
        removeFile('optics_'+observationID+'.pars')
        removeFile('catlist_'+observationID+'.pars')

        if self.grid == 'no':
            for p in jobs:
                p.join()
        elif self.grid == 'condor':
             condor.submitDag(self)
        elif self.grid == 'diagrid':
             diagrid.submitDax(self)
        os.chdir(self.phosimDir)
        return rImageArg

     ## Generic methods for handling execution environment
     def initExecutionEnvironment(self):
        if self.execEnvironmentInitialized:
            return
        if self.grid == 'condor':
            self.initCondorEnvironment()
        elif self.grid == 'diagrid':
            self.initDiagridEnvironment()
        elif self.grid == 'cluster':
            self.initClusterEnvironment()
        self.execEnvironmentInitialized = True

     ## general method to delete files at end
     def cleanup(self,keep_screens,rImageArg):
        if self.grid in ['no', 'cluster']:
            os.chdir(self.workDir)
            removeFile('objectcatalog_'+self.observationID+'.pars')
            removeFile('tracking_'+self.observationID+'.pars')
            removeFile('fea_'+self.observationID+'_*.txt')
            if not keep_screens:
                 removeFile('airglowscreen_'+self.observationID+'.fits.gz')
                 for f in glob.glob('atmospherescreen_'+self.observationID+'_*') :
                      removeFile(f)
                 for f in glob.glob('cloudscreen_'+self.observationID+'_*') :
                      removeFile(f)
            else:
                 f='atmosphere_'+self.observationID+'.pars'
                 shutil.move(f,self.outputDir+'/'+f)
                 f='airglowscreen_'+self.observationID+'.fits.gz'
                 shutil.move(f,self.outputDir+'/'+f)
                 for f in glob.glob('atmospherescreen_'+self.observationID+'_*') :
                      shutil.move(f,self.outputDir+'/'+f)
                 for f in glob.glob('cloudscreen_'+self.observationID+'_*') :
                      shutil.move(f,self.outputDir+'/'+f)
            if self.eventfile==1:
                 for f in glob.glob('*_r_'+self.observationID+'_*'):
                      shutil.move(f,self.outputDir+'/'+f)
            if self.throughputfile==1:
                 for f in glob.glob('throughput_*'+self.observationID+'_*') :
                      shutil.move(f,self.outputDir+'/'+f)
            if self.centroidfile==1:
                 for f in glob.glob('centroid_*'+self.observationID+'_*') :
                      shutil.move(f,self.outputDir+'/'+f)
            for f in glob.glob('opd*'+self.observationID+'_*') :
                 shutil.move(f,self.outputDir+'/'+f)
            os.chdir(self.phosimDir)

            if rImageArg != '':
                 commandName = '${PYTHON} '+os.path.join(self.phosimDir, 'tools/phosim_visualizer/phosim_visualizer.py')+' '+self.instrDir+'/ '+rImageArg
                 if subprocess.call(commandName, shell=True) != 0:
                      sys.exit(1)

     ## Condor method to setup directories
     def initCondorEnvironment(self):
        sys.path.append(self.phosimDir+'/condor')
        global condor
        import condor
        condor.initEnvironment(self)

     def initDiagridEnvironment(self):
        sys.path.append(self.phosimDir+'/diagrid')
        global diagrid
        import diagrid
        diagrid.initEnvironment(self)

     ## Cluster methods
     def initClusterEnvironment(self):
        pass

##main function
def main():

     phosimDir=os.path.split(os.path.abspath(__file__))[0]
     defaultOutputDir=os.getenv('PHOSIM_OUTPUT_DIR', phosimDir+'/output')
     defaultWorkDir=os.getenv('PHOSIM_WORK_DIR', phosimDir+'/work')
     defaultBinDir=os.getenv('PHOSIM_BIN_DIR', phosimDir+'/bin')
     defaultDataDir=os.getenv('PHOSIM_DATA_DIR', phosimDir+'/data')
     defaultSedDir=os.getenv('PHOSIM_SED_DIR', phosimDir+'/data/SEDs')
     defaultImageDir=os.getenv('PHOSIM_IMAGE_DIR', phosimDir+'/data/images')

     cpuCount = 1
     try:
        cpuCount = multiprocessing.cpu_count()
     except NotImplementedError:
        pass
     
     parser = optparse.OptionParser(usage='%prog instance_catalog [<arg1> <arg2> ...]')
     parser.add_option('-c','--command',dest="extraCommands",default="none",
             help='command file to modify the default physics')
     parser.add_option('-p','--proc',dest="numproc",default=1,type="int",
             help='number of processors')
     parser.add_option('-t','--thread',dest="numthread",default=cpuCount,type="int",
             help='number of threads')
     parser.add_option('-o','--output',dest="outputDir",default=defaultOutputDir,
             help='output image directory')
     parser.add_option('-w','--work',dest="workDir",default=defaultWorkDir,
             help='temporary work directory')
     parser.add_option('-b','--bin',dest="binDir",default=defaultBinDir,
             help='binary file directory')
     parser.add_option('-d','--data',dest="dataDir",default=defaultDataDir,
             help='data directory')
     parser.add_option('--sed',dest="sedDir",default=defaultSedDir,
             help='SED file directory')
     parser.add_option('--image',dest="imageDir",default=defaultImageDir,
             help='truth image directory')
     parser.add_option('-s','--sensor',dest="sensor",default="all",
             help='sensor chip specification (e.g., all, R22_S11, "R22_S11|R22_S12")')
     parser.add_option('-i','--instrument',dest="instrument",default="lsst",
             help='instrument site directory')
     parser.add_option('-g','--grid',dest="grid",default="no",
             help='execute remotely (no, condor, cluster, diagrid)')
     parser.add_option('-u','--universe',dest="universe",default="standard",
             help='condor universe (standard, vanilla)')
     parser.add_option('-e','--e2adc',dest="e2adc",default=1,type="int",
             help='whether to generate amplifier images (1 = true, 0 = false)')
     parser.add_option('--keepscreens',dest="keepscreens",default=0,type="int",
             help='whether to keep atmospheric phase screens (0 = false, 1 = true)')
     parser.add_option('--checkpoint',dest="checkpoint",default=12,type="int",
             help='number of checkpoints (condor only)')
     parser.add_option('--ds9',dest="ds9",action="store_true",
             help='whether to launch ds9 upon completion')
     parser.add_option('--visualize',dest="visualize",action="store_true",
             help='whether to launch phosim visualizer upon completion')
     parser.add_option('-v','--version',action='store_false',help='prints the version')

     if len(sys.argv)<2:
        usage()
        sys.exit()

     if sys.argv[1] in ('-h', '--help'):
        usage()
        sys.exit()

     if sys.argv[1] in ('-v', '--version'):
        print(' ')
        print('Photon Simulator (PhoSim)')
        os.system('more '+defaultBinDir+'/version')
        print(' ')
        sys.exit()

     opt, remainder = parser.parse_args(sys.argv[1:]) #parse_args returns a pair of values
     instanceCatalog=remainder[0]

     checkPaths(opt, phosimDir)

     grid_opts = {'numproc': opt.numproc}
     if opt.grid == 'condor':
          grid_opts = {'universe': opt.universe, 'checkpoint': opt.checkpoint}
     elif opt.grid == 'diagrid':
          grid_opts = {'checkpoint': opt.checkpoint}
     elif opt.grid == 'cluster':
          grid_opts = {'script_writer': jobChip}

     #the entire phosim workflow follows:
     focalplane = PhosimFocalplane(phosimDir, opt, grid_opts, bool(opt.visualize))
     focalplane.loadInstanceCatalog(instanceCatalog, opt.extraCommands)
     os.chdir(opt.workDir)
     focalplane.writeInputParamsAndCatalogs()
     focalplane.generateAtmosphere()
     focalplane.generateInstrumentConfig()
     focalplane.generateSimulatedCatalogs()
     focalplane.trimObjects(opt.sensor)
     rImageArg = focalplane.scheduleRaytrace(opt.instrument, bool(opt.e2adc),bool(opt.keepscreens),bool(opt.ds9),bool(opt.visualize))
     focalplane.cleanup(bool(opt.keepscreens),rImageArg)

if __name__ == "__main__":
    main()
