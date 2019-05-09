#!/usr/bin/env python
##
## @package phosim.py
## @file cluster_submit.py
## @brief script to submit phosim jobs to a cluster
##
## @brief Created by:
## @author Glenn Sembroski (Purdue)
##
## @brief Modified by:
##
## @warning This code is not fully validated
## and not ready for full release.  Please
## treat results with caution.
##
## The intention of this script is to use the dagMan and  .submit files
## created in the phosim "work" directory by the running of phosim with the 
##  -g condor option, to generate and submit the slurm (or pbs) files to run 
## the trim, raytrace and e2adc jobs on a cluster. And then to submit them.
## Note: Particular clusters required specific modifications to the submission
## files. (Provisions and examples of this exist in functions
## setupForHost() and getPid(0 functions).
#

import os
import sys, optparse
import subprocess
from collections import OrderedDict   #For getJobs

phosimVersion = 'dev'


## print the usage
def usage():
     global phosimVersion
     script=os.path.abspath(__file__)
     os.system(script+' x --help')
     print('Phosim Version: ' + phosimVersion)
     return

##For debuggign we fake a sumssion return to get a PID
def generateDebugSubPBSLine(jobCount):
    host = os.getenv('HOSTNAME')
    j= str(jobCount)
    if host[:6] == "hammer":
        subPBSLine = (j + '.hammer-adm.rcac.purdue.edu')
    else:
        subPBSLine = ('Submitted batch job ' + j )
    return subPBSLine

## Seach the dagMan file for the dependancy lines. They start with the word
# "PARENT". Create and return a dict of Child:Parent for use later is
# deriveing dependancy commands for the submission file
def getDependancies(dagManFileFull):
     #search the dagMan file for the "PARENT" lines. On these lines are the
     #job dependancies

     #Get all the lines
     dFile = open(dagManFileFull,"r")
     dagMan = dFile.readlines()
     dFile.close()

     dependancies = {}
     for line in dagMan:
          lstr1 = line.split()
          #see if this line is a depancncy line
          if lstr1[0] == "PARENT" and lstr1[2] == "CHILD":
               #extract things from this line
               parentJobName = lstr1[1]
               childJobName = lstr1[3]
               dependancies[childJobName] = parentJobName
     return dependancies

##getJobs returns an OrderedDict with the jobName as key and submission file name
# as the value
def getJobs(dagManFileFull,jobType,requiredExtention=None):

    #Go through the dagman file and find all "JOB" lines with jobType
    #and, if specfied, end in requiredExtention.

    #Get all the lines (I know this duplicates the reading of the file we did
    #in getDependancies but I wanted to keep thing independant and the code
    #obvious
     
    dFile = open(dagManFileFull,"r")
    dagMan = dFile.readlines()
    dFile.close()
     
    #Search for "JOB" and jobType lines (with requiredExtention if specified)
    #This is to ingnore in the condor --checkpointcount =1 (or some
    #value !=0) generated dagman file job lines. The users should really
    #specifiy checkpointing when running phosim but we handle it (and ignore
    #it if he does. Checkpointing is slected by using the checkpoint option
    #when running this script.
    #If we are going to have internal checkpointing we assume the order of the
    #jobs presented in the dagman files is the order they should be submitted
    #(so we can propogate the  dependencies properly) Thus we use OrderedDict
    #here.
    
    jobs = OrderedDict()
    jobTypeLen=len(jobType)
    for line in dagMan:
        lstr1 = line.split()
        #see if this is a JOB and jobType and requiredExtention line
        if len(lstr1) > 1 :                       #Ignore blank lines
            if lstr1[0] == "JOB":                 #A JOB line
                jobName = lstr1[1]
                if jobName[:jobTypeLen] == jobType :#A JobType Line
                    if requiredExtention is None :  #Test for requiredExtention
                        goodExtention=True
                    else:
                        extenLen = len(requiredExtention)
                        goodExtention = (jobName[-extenLen:] ==
                                                              requiredExtention)
                    if goodExtention : 
                        submitFileName = lstr1[2]
                        jobs[jobName] =  submitFileName
    return jobs

##Read througth the submision file for this job and extract the values into
# a dict.
def getSubmissionParams(submitFileName):

    sFile = open(submitFileName,"r")
    jobSubmitLines = sFile.readlines()
    sFile.close()
               
    jobSubDict = {}
    for line in jobSubmitLines:
        submitLine=line.split()
        if submitLine[0] == "executable":
            jobSubDict['executable'] = submitLine[2]
        elif submitLine[0] == "initialdir":
            jobSubDict['initialDir'] = submitLine[2]
        elif submitLine[0] == "Input":
            jobSubDict['input'] = submitLine[2]
        elif submitLine[0] == "Output":
            jobSubDict['output'] = submitLine[2]
        elif submitLine[0] == "Log":
            jobSubDict['log'] = submitLine[2]
        elif submitLine[0] == "Error":
            jobSubDict['error'] =  submitLine[2]
        elif submitLine[0] == "transfer_input_files":
            jobTransferFiles = " "
            jobTransferFiles = jobTransferFiles.join(submitLine[2:])
            jobSubDict['transferFiles'] = jobTransferFiles

    return jobSubDict

## Setup a string of pbs/slrumn values;
def setupForHost(knl,haswellShared,haswellRegular):
    host = os.getenv('HOSTNAME')
    #Setup the commands for the cluster we are using:
    #Setup the nodes and cpu per node. We will change this a little when we 
    #add the possibility of threads. 2 possible options either seperate nodes 
    #and ppn or combined(NERSC)
    #   
    
    #For NERSC we may want to list the disk we are using.
    #Returend is: base of PBS/slrumn submit command,Dependency command,
    #string of all other commands, string of maximum thread count.
    if haswellShared :
         print('Haswell:Shared: True')
    elif haswellRegular:
         print('Haswell:Regular: True')

    if host[:4] != 'cori':
          knl = False
          haswellRegular = False
          haswellShared = False
          
    if host[:6] == 'edison':
        pbsSetupList = ('#SBATCH -L SCRATCH,project' + "\n" +
                        '#SBATCH -p shared' + "\n" +
                        '#SBATCH -t 10:00:00' + "\n" +
                        '#SBATCH -N 1' + "\n" +
                        '#SBATCH --mem=5GB' + "\n" +
                        '#SBATCH -A m1727'+ "\n")

        return {'HOST':          'edison',
                'SUBMITCMD':     'sbatch',
                'SUBMITJOBID':   'SLURM_JOB_ID',
                'DEPENDCMD':     '#SBATCH -d afterok:', 
                'INITIALLIST':   pbsSetupList,
	        'MAXTHREADS':    '24',
                'THREADCMD':     '#SBATCH -n ',
                'DMTCPSRC':      '/global/common/edison/contrib/lsst/dmtcp',
                'DMTCPVERSION':  'dmtcp-2.5.0',
                'DMTCPINTERVAL': '3600'}

    elif host[:4] == 'cori' and  haswellShared :
        pbsSetupList = ('#SBATCH -L SCRATCH,project' + "\n" +
                        '#SBATCH -p shared' + "\n" +
                        '#SBATCH -t 40:00:00' + "\n" +
                        '#SBATCH -N 1' + "\n" +
                        '#SBATCH -C haswell' + "\n" +
                        '#SBATCH -A m1727'+ "\n")
 
        return {'HOST':          'cori',
                'SUBMITCMD':     'sbatch',
                'SUBMITJOBID':   'SLURM_JOB_ID',
                'DEPENDCMD':     '#SBATCH -d afterok:', 
                'INITIALLIST':   pbsSetupList,
                'MAXTHREADS':    "32",
                'THREADCMD':     '#SBATCH -n ',
                'DMTCPSRC':      '/global/common/cori/contrib/lsst/dmtcp',
                'DMTCPVERSION':  'dmtcp-2.5.0',
                'DMTCPINTERVAL': '3600'}   #In seconds

    elif host[:4] == 'cori' and  haswellRegular :
        pbsSetupList = ('#SBATCH -L SCRATCH,project' + "\n" +
                        '#SBATCH -p regular' + "\n" +
                        '#SBATCH -t 40:00:00' + "\n" +
                        '#SBATCH -N 1' + "\n" +
                        '#SBATCH -C haswell' + "\n" +
                        '#SBATCH -A m1727'+ "\n")
 
        return {'HOST':          'cori',
                'SUBMITCMD':     'sbatch',
                'SUBMITJOBID':   'SLURM_JOB_ID',
                'DEPENDCMD':     '#SBATCH -d afterok:', 
                'INITIALLIST':   pbsSetupList,
                'THREADCMD':     '#SBATCH -n ',
                'MAXTHREADS':    "60",
                'DMTCPSRC':      '/global/common/cori/contrib/lsst/dmtcp',
                'DMTCPVERSION':  'dmtcp-2.5.0',
                'DMTCPINTERVAL': '3600'}   #In seconds
   
    elif host[:4] == 'cori' and  knl :
        pbsSetupList = ('#SBATCH -L SCRATCH,project' + "\n" +
                        '#SBATCH -p regular' + "\n" +
                        '#SBATCH -t 10:00:00' + "\n" +
                        '#SBATCH -N 1' + "\n" +
                        '#SBATCH -C knl,quad,cache' + "\n" +
                        '#SBATCH -S 4       # use core specialization' + '\n'
                        '#SBATCH -A m1727'+ "\n")
 
        return {'HOST':          'cori',
                'SUBMITCMD':     'sbatch',
                'SUBMITJOBID':   'SLURM_JOB_ID',
                'DEPENDCMD':     '#SBATCH -d afterok:', 
                'INITIALLIST':   pbsSetupList,
                'MAXTHREADS':    "268",
                'THREADCMD':     'SLURM Thread command undefined for KNL!' + "\n",
                'DMTCPSRC':      '/global/common/cori/contrib/lsst/dmtcp',
                'DMTCPVERSION':  'dmtcp-2.5.0',
                'DMTCPINTERVAL': '3600'}   #In seconds

    elif host[:6] == 'hammer':
        #pbsSetupList = ('#PBS -q physics' + "\n" +
        #                '#PBS -l walltime=10:00:00' + "\n" +
        pbsSetupList = ('#PBS -q standby' + "\n" +
                        '#PBS -l walltime=04:00:00' + "\n" +
                        '#PBS -l mem=8GB' + "\n" +
                        '#PBS -l naccesspolicy=singleuser'+ "\n")

        return {'HOST':          'hammer',
                'SUBMITCMD':     'qsub -V ',
                'SUBMITJOBID':   'PBS_JOBID',
                'DEPENDCMD':     '#PBS -W depend=afterok:', 
                'INITIALLIST':   pbsSetupList,'MAXTHREADS':"20",
                'THREADCMD':     '#PBS -l nodes=1:ppn=',
                'DMTCPSRC':      '/depot/lsst/apps/dmtcp',
                'DMTCPVERSION':  'dmtcp-2.5.0',
                'DMTCPINTERVAL': '3600'}

    elif host[:8] == 'halstead':
        #pbsSetupList = ('#PBS -q physics' + "\n" +
        #                '#PBS -l walltime=10:00:00' + "\n" +
        pbsSetupList = ('#PBS -q standby' + "\n" +
                        '#PBS -l walltime=04:00:00' + "\n" +
                        '#PBS -l mem=8GB' + "\n" +
                        '#PBS -l naccesspolicy=singleuser'+ "\n")

        return {'HOST':          'halstead',
                'SUBMITCMD':     'qsub -V ',
                'SUBMITJOBID':   'PBS_JOBID',
                'DEPENDCMD':     '#PBS -W depend=afterok:', 
                'INITIALLIST':   pbsSetupList,'MAXTHREADS':"20",
                'THREADCMD':     '#PBS -l nodes=1:ppn=',
                'DMTCPSRC':      '/depot/lsst/apps/dmtcp',
                'DMTCPVERSION':  'dmtcp-2.5.0',
                'DMTCPINTERVAL': '3600'}

    else:
        print('Unknow host: ' + host )
        print('Known hosts are: NERSC:edison, cori (haswell or KNL), ' +
              'Purdue: hammer, halstead')
        sys.exit(1)

##Search Input pars file for "thread" command and get the number of requested
#threads (as a string). If not there return requested number of "1"
def getNumRequestedThreads(inputFileName):
    iFile = open(inputFileName,'r')
    iLines = iFile.readlines()
    iFile.close()
    thread=1
    for line in iLines:
        lstr1 = line.split()
        #see if this line is a "thread" line
        if len(lstr1) == 2 :
            if lstr1[0] == "thread" :
                #extract thread count  from this line
                threadCount = lstr1[1]
                break

    return threadCount

##Search Input pars file for "checkpoint" command. Return True if
#checkpointing requested 
def isPhosimChckptRequested(inputFileName):
    iFile = open(inputFileName,'r')
    iLines = iFile.readlines()
    iFile.close()
    checkpointRequested = False 
    for line in iLines:
        lstr1 = line.split()
        #see if this line is a "checkpointtotal" line
        if len(lstr1) == 2 :
            if lstr1[0] == "checkpointtotal" :
                #extract checkpoint request status (0 = False, 1 = True)
                cpStatus = lstr1[1]
                if cpStatus != "0" :
                    checkpointRequested = True
                break

    return checkpointRequested

##Copy the raytrace .pars file to a new file disabling the internal raytrace
#checkpointing options. Set both checkpointtotal and checkpointcount lines to 0
#Returns name of new .pars file (which we set here to be "nocp" + parsFile)
def createNonCheckpointingParsFile(inputDir, parsFile):
    #if checkpointtotal not 0 add the sed commands to set the chkpointtotal and
    #checkpointcount to 0.

    parsFileFull = inputDir + '/' + parsFile
    internalCheckPointReq = isPhosimChckptRequested(parsFileFull)
    if not internalCheckPointReq:
        return parsFile
    
    inputFile = open(parsFileFull, 'r')

    nocpParsFile = "nocp" + parsFile
    nocpParsFileFull =  inputDir + '/' + nocpParsFile
    outputFile = open(nocpParsFileFull, 'w')

    for line in inputFile:
        lstr1=line.split()
        #Note: we remove blank lines
        if len(lstr1) >0 :
            if lstr1[0] == 'checkpointtotal' :
                outputFile.write('checkpointtotal 0' + '\n')
            elif lstr1[0] == 'checkpointcount' :
                outputFile.write('checkpointcount 0' + '\n')
            else:
                outputFile.write(line)

    inputFile.close()
    outputFile.close()
    return nocpParsFile 

##Add bash commands and defintions to submit file to setup for running
#with DMTCP. Create (or make sure it exists) the dmtcp and dmtcp/tmp* 
#directories off of "work"
#This could also be the place where we check that dmtcp exists on this cluster:
#assume it does.
def addDMTCPSetupCommands(pfile, submitPBSList, jobName) :
    #Create (if they done't exist the work/dmtcp/tmp* directories
    #assume we have already cd to work direcdtory in submit file.
    dmtcpTmpDirName = ('tmp' + jobName[9:-2] + '_${' +
                submitPBSList["SUBMITJOBID"] + '}' ) 
    pfile.write("mkdir -vp dmtcp/" + dmtcpTmpDirName + '\n')

    #setup the export definition commands for dmtcp
    pfile.write("lcl=$PWD" + '\n')
    pfile.write("export DMTCP_TMPDIR=$lcl/dmtcp/" + dmtcpTmpDirName  + '\n')
    pfile.write("export DMTCP_CHECKPOINT_DIR=$DMTCP_TMPDIR" + '\n')
    ## activates auto ckpt (in seconds)
    pfile.write("export DMTCP_CHECKPOINT_INTERVAL=" +
                                        submitPBSList["DMTCPINTERVAL"] + '\n')

    pfile.write("export DMTCP_ROOT=" + submitPBSList["DMTCPSRC"] + '/'
                +submitPBSList["DMTCPVERSION"] + '\n')
    pfile.write("export PATH=$DMTCP_ROOT/bin:$PATH" + '\n')
    pfile.write("export MANPATH=$DMTCP_ROOT/share/man:$MANPATH" + '\n')

    return

## Get the cluster assiged PID for this submitted job.
#This usually only called for the trim jobns and used when setting up 
#dependencies for the raytrace /e2adc jobs.
def getPid(jobPBS,subPBSLine):
    #print ('Submission of ' + jobPBS  + ' returned:' + 
    #       subPBSLine)
    host = os.getenv('HOSTNAME')
    if host[:6] == 'edison':
        #return from the submission will be something like:
        # "Submitted batch job 2162058". So we want the 4th word

        JobPid = subPBSLine.split()[3]

    if host[:4] == 'cori':
        #return from the submission will be something like:
        # "Submitted batch job 2162058". So we want the 4th word

        JobPid = subPBSLine.split()[3]

    if host[:6] == 'hammer':
        #return from the submission will be something like:
        # "1029466.hammer-adm.rcac.purdue.edu". So we want the string before 
        #the "."

        JobPid = subPBSLine.split(".")[0]

    if host[:8] == 'halstead':
        #return from the submission will be something like:
        # "1029466.halstead-adm.rcac.purdue.edu". So we want the string before 
        #the "."

        JobPid = subPBSLine.split(".")[0]
   
    return JobPid
##On Cori (knl or haswellShared or haswellRegular will be true to get here) add
#the module load and module list commands
def addModuleLoadCommands(pfile,knl, haswellShared, haswellRegular):
     global phosimVersion
     if knl :
          pfile.write('module load phosim/' + phosimVersion + 'KNL' + '\n')
          #pfile.write('module load phosim/' + phosimVersion + 'KNLFast1' + '\n')
     elif haswellRegular or haswellShared :
          pfile.write('module load phosim/' + phosimVersion + 'Haswell' + '\n')
         
     pfile.write('module list' + '\n')
     return

## Add the sbatch (or pbs) cluster commands to the submission file.
def addClusterCommands(pfile, submitPBSList, threadCMD, numThreads,knl,
                       haswellRegular):
         #Print out all the PBS/slrumn commands in the List (except for the 
    #DEPEND cmd. Thats later.
    #The submitPBSList and threadCMD strings come from the call to 
    #setupForHost 
    #Setup start of pbs submission file
    
    pfile.write("#!/bin/bash -l\n")
    pfile.write(submitPBSList)
    if ( not knl ) and (not haswellRegular) :
         pfile.write(threadCMD + numThreads + "\n")
    return

## Determine if this is "CHILD" job. That is if it depends on a parent job
# running first. If so add a depenacny line to the submission file.
def addDependancyCommand(pfile, jobName, dependancies, jobIDs, 
                         submitDependCmd):
    # First check to see that this jobName is a child to someones parent
    # Then check that the parent has already been run and we have a pid for it
    #
    # This could be more general by not assuming the above ordering by
    # searching for jobs with no child tasks, running them first and so on
    # and so on. That would complicate thing highly but its doable
    if dependancies.get(jobName) is not None:
        parentJobName=dependancies[jobName]           #Parent trim job name

        #print('parentJobName,jobName: '+ parentJobName + ' ' + jobName + '\n')
        #At this point things get complicated.
        #Now make sure that this parent has run and if so get its pid
        if jobIDs.get(parentJobName) is not None:

            pid =  jobIDs[parentJobName]
            # Add the dependancy command
            pfile.write(submitDependCmd + pid + '\n')
        else:
            print('No pid avaiable for dependancy ' + parentJobName +
                  ' needed by ' +jobName )
            sys.exit(1)
    return

## Add to the pbs file the lines to run this task. (Originally from the dagMan 
# file.
def addTaskRunLine(pfile, jobSubDict, taskName, parsFile, dmtcpChckPointReq,
                   appendToLog = False):
    #run the job with its arguments. This will be multiline with 
    #continuations.
    #The log file will need to manipulated  a bit to start with the taskname
    path,logFileName = os.path.split(jobSubDict["log"])
    logFileFull = (jobSubDict["initialDir"] + "/" + path + '/' + taskName +
                   logFileName[3:])
    # For dmtcp checkpointing: Preface the task run command with the
    # dmtcp_launch command with its options
    if  dmtcpChckPointReq :
        pfile.write("$DMTCP_ROOT/bin/dmtcp_launch --new-coordinator " +
                     "--port-file  $DMTCP_TMPDIR/port.txt " + " \\"  + '\n')

    pfile.write("time " + jobSubDict['executable'] + " \\"  + '\n')
    pfile.write("< " + parsFile + " \\"  + '\n')
    if appendToLog:
         pfile.write(">>" + logFileFull + '\n')
    else:
         pfile.write(">" + logFileFull + '\n')

    return

##Add a line to the .pbs script file to remove a file. 
def addRemoveFileLine(pfile,fileName): 
    pfile.write('rm ' + fileName + '\n')
    return

#Add the commands to the submission job script that runs the job.
# Not sure about the transfer stuff for now. Wait and see whats missing
def addJobCommands(pfile,jobSubDict,taskName,jobName,dmtcpChckPointReq,
                   submitPBSList,parsFileName):
     
    # #############################################################
    # Note: We do a lot of the manipulating of the various input and files
    # (for raytrace this is:checkpointing, trim catalog etc) dynamically in
    # the submission file, which gets executed on the compute node at run time.
    # #############################################################
    
    #Add the bash commands to the submission file that sets us up for running
    #this job.
    
    #On First set to write start time and to cd to the initial directory
    pfile.write('date' + '\n')
    pfile.write('cd ' + jobSubDict['initialDir'] + '\n')

    #Trim jobs:No checkpointing, new log file.
    if taskName != 'raytrace':    
        addTaskRunLine(pfile, jobSubDict, taskName, parsFileName, False, False)
        return

    elif taskName == 'raytrace':
        #Now add in the export bash commands to set the dmtcp control
        #values
        if dmtcpChckPointReq :
            addDMTCPSetupCommands(pfile, submitPBSList, jobName)    

        #Because the raytrace/e2adc jobs can be run multiple times we need to
        #make a temperary raytrace pars file each time (so we only get one 
	#copy of the trim catalog in it.
        #So:
	#1: Add a line in the submission file that copies (makes a new one) 
	#   the raytrace*.pars file to a new file with a temp name (*.tmp).
	#2: Add a line in the submission file to append the trimcatalog (which 
	#   the trim jobs make) to the end of the temp raytrace.pars file
        #Trim catalog and raytrace.pars files will be in local directory(work)

	tempParsFileName = parsFileName + '.tmp'
	pfile.write( 'cp ' +  parsFileName + ' ' + tempParsFileName + '\n')

        trimCatalogFileName = 'trimcatalog_' + jobName[9:-7] + '.pars'
        pfile.write( 'cat ' + trimCatalogFileName + ' >> ' +
	             tempParsFileName + '\n')

        #Determine if this is the first raytrace of this chip where a new log
        #file will be needed. Internal phosim checkpoint runs will have
        #additional dependent raytrace jobs whose log ouputs we will want to
        #append to the intial raytrace log file.
        if jobName[-2:] == "_0":
            appendToLog = False
        else:
            appendToLog = True
            
        #I may need to transfer or make links to some files here. Not sure
        addTaskRunLine(pfile, jobSubDict, taskName, tempParsFileName,
                       dmtcpChckPointReq,appendToLog)
    return


## removeFile deletes files (if they do not exist, it will catch the OSError 
# exception and silently proceed.)
def removeFile(filename):
    try:
         os.remove(filename)
    except OSError:
         pass

## Add a line to a .pbs file to mv a file to a directory
def addMoveFileLine(pfile,fileName,destDir):
    pfile.write( 'mv ' + fileName + ' ' + destDir + '/\n')
    return

## Seach the dagMan file for the jobs. Create the slurm ( or .pbs
# depending on cluster)
# submit the jobs and save JobNames and cluster Job_ID in a list which is
# returned. This list will be used for dependancies
def createAndSubmitJobs(opt,dagManFileFull):
     
    #Check that both knl and haswell partions haven't both been selected
    global phosimVersion
    knl =  opt.coriKNLPartition
    phosimVersion = opt.PhosimVersion
    haswellShared = opt.coriHaswellSharedPartition
    haswellRegular = opt.coriHaswellRegularPartition
    if knl and ( haswellShared or haswellRegular):
         print("False-- Both KNL and HASWELL paritions selected!")
         sys.exit(1)
    elif ( not knl)  and (not haswellRegular) and (not haswellShared) :
         haswellShared = True       #Default to haswellShared partition
         
 
    #Setup the cluster specific stuff.
    submitPBSList = setupForHost(knl, haswellShared, haswellRegular)
    
    #Go through the dagman file and make a dict of PARENT:CHILD jobname
    #dependancies.
    dependancies = getDependancies(dagManFileFull)


               
    #Search the dagMan file for the jobs. We use a suggestion from En-hsin ping
    #here: We assume that the order of all jobs in the dagman file is in the
    #order of dependencies when it comes to the phosim internal checkpointing.
    #So we use the OrderedDict class.
    #getJobs returns a OrderedDict with the jobName as key and submission file 
    #name as the value.  These dictionaries just relate a job name with
    #its .submit file name which has all the input pars filename in it among
    #other useful things     
    trimJobDict = getJobs(dagManFileFull,"trim_")
    raytraceJobDict = getJobs(dagManFileFull,"raytrace_")
    raytraceJobDict_0 = getJobs(dagManFileFull,"raytrace_","_0")

    e2adcJobDict =  getJobs(dagManFileFull,"e2adc_")

    # ####################################################################
    #Get flags for various modes of running and test that the combinations make
    #sense
    disableTrimSubmission = opt.disableTrim
    disableRaytraceSubmission = opt.disableRaytrace
    disableE2adcSubmission = opt.disableE2adc

     #Make sure they make sense(see comments below)
    if disableTrimSubmission and  disableRaytraceSubmission : 
         print("WHAT? Fatal--No use running this script if its not going to "
               "do anything!" + '\n')
         sys.exit(1)

    #The disableE2adcSubmission flag only has validity(makes sens e to use) if
    #we are making up new raytarce jobs. That can only happen if we are
    #submitting trim jobs(we need the pid's)
    if disableTrimSubmission and disableE2adcSubmission:
         print("WARNING-   --E2adc-Disable option invalid when used with "
               "--Trim-Disable option!" + '\n')
 

    # ######################################################################
    #jobNamePID is a dict used for generating raytrace*_n job PID
    #dependancies (no matter which if any checkpointing was requested).
    #jobNamePID links the trim (or raytrace)  JobName and submitted job Pid's
    #(process id). This is created from the jobs submissions only and
    #used only for setting the dependiencies of the raytrace*_n jobs
    #The trim jobs will have no dependencies and each 
    #raytrace*_0 job will only depend on a particular trim job (many raytrace 
    #jobs (9*NumExposures acturally, a raft) will depend on a particulat trim
    #job finishing). Each e2adc job depends on a particular raytrace_n job
    #(which raytrace job depends
    #on if we are doing internal phosim checkpointuing or not. We implictly
    #enforce that dependency by including the running of the e2adc task in its 
    #Parent raytrace_n pbs job file.
    jobNamePID={}

    #Trim jobs first(Order follows order in dagman file since trimJobDict is
    #an OrderedDict
    jobCount = 0

    #Skip creating and submitting trim jobs if options tells us to.
    if not disableTrimSubmission:
        for jobName,trimSubmitFileName in trimJobDict.items():	
            #Get a dict of the important lines from the submission file
            jobSubDict = getSubmissionParams( trimSubmitFileName)
                  
            #We have everything we need to submit this task.
            #time to build the .pbs file
            # ##############################
            #Note For the purposes of this script the .pbs abreviation is used 
            #for BOTH pbs AND  slurm based clusters)
            # ##############################

            jobPBS = opt.workDir + "/" + jobName + '.pbs'
            pfile=open(jobPBS,'w') 
            numTrimThreads = "1"     #Note this is a string!  
            addClusterCommands(pfile, submitPBSList['INITIALLIST'],
                               submitPBSList['THREADCMD'], numTrimThreads,
                               knl, haswellRegular)
            addModuleLoadCommands(pfile,knl, haswellShared, haswellRegular)
            #No dmtcp checkpointing for trim jobs(though we could if needed)
            #Add the trim command.
            #The 'log' just retains trim log file name as listed in the .submit 
            #file.  That name alteady has 'trim' in it
            trimParsFileName = jobSubDict['input']
            addJobCommands(pfile, jobSubDict, 'log', jobName, False,
                           submitPBSList, trimParsFileName)
  
            #cleanup(remove .pars file from "work/" Disabled for now.
            #addRemoveFileLine(pfile, jobSubDict['input'])
                    
            # The next line will probably work but lets wait a bit before
            # we enable it. Disabled for now.
	        #addRemoveFileLine(pfile,jobPBS)

            #And thats all for trim
            pfile.close()
          
            #submit the job
            #Its important here to get the Job_ID (PID) that gets returned to
            #standard output for use later in setting the raytrace job
            #dependency.
          
            submitCmd = (submitPBSList['SUBMITCMD'] + 
                         ' -e ' + jobSubDict['initialDir'] + '/errors/' +
                         jobName + '.pbs.err ' + 
                         ' -o ' + jobSubDict['initialDir'] + '/logs/' +
                         jobName + '.pbs.log ' + jobPBS)

            #print ('Submit cmd: ' + submitCmd) 

            #This is it! RUN iT!!!
            # #################################################
        
            subPBSLine = subprocess.check_output(submitCmd, shell=True)
        
            #Debug lines follows
            #jobCount = jobCount +1
            #subPBSLine = ('Submitted batch job ' + str(jobCount) )

            print(jobName + ': ' + subPBSLine)

            #Pick up the PID for dependencies later
            # We should probably check the return code here. Later!
            JobPid = getPid(jobPBS,subPBSLine)
            # #################################################
               


            #save this in a dict with the trim job name as a key. Used
            # when setting up raytrace dependencies
            jobNamePID[jobName]=JobPid

            #print('jobName,JobPid: ' + jobName + ' ' + JobPid + '\n')

            #cleanup
            #removeFile(trimSubmitFileName) Disabled for now.

        #If enabled the Trim jobs are all now all running and we have the PID
        # values for the dependency commands.

    # ################################################################
    #Build and submit the combined raytrace jobs.
    #To make our life easy we have built (above) a dict with the raytrace*_0
    #JOBS. 
    #We now determine the number of threads from the first raytrace job's
    #(probably a raytrace*_0) .pars file.
    #Get jobName (Key) and submit file name
    jobName,raytraceSubmitFileName = raytraceJobDict.items()[0]
    jobSubDict = getSubmissionParams(raytraceSubmitFileName)
    inputDirFull = (jobSubDict["initialDir"] + "/" +  jobSubDict["input"])

    ##############THREADS########
    #get the requested threads(if asked for or max if exceeded)
    requestedThreads = getNumRequestedThreads( inputDirFull )
    if int(requestedThreads) > int(submitPBSList['MAXTHREADS']) :
         requestedThreads=submitPBSList['MAXTHREADS']

    #########INTERNAL CHECKPOINTING##################
    #Determine type of checkpointing, if any
    #Two possible types : phosim internal and dmtcp checkpointing (kind of
    #external) From command line:
    phosimChckPointReq = opt.phosimChckpntReq      
    dmtcpChckPointReq = opt.dmtcpChckpntReq

    #dmtcp Checkpointing will take priority at this time (If you want both we
    #can probably do that too)
    if dmtcpChckPointReq :
        #Check we are on a node where dmtcp works
        if (submitPBSList['HOST'] == 'cori' or
            submitPBSList['HOST'] == 'edison' or
            submitPBSList['HOST'] == 'hammer' ) :
             phosimChckPointReq = False     #Force no internal checkpointing
        else:
             print('--dmtcp-checkpointing ( -d ) option using dmtcp, only ' +
                   'available for clusters cori,edison,hammer ' +
                   'at present')
             print('No dmtcp checkpointing implimented')
             dmtcpChckPointReq = False
             
    #PhoimCheckPointReq only true at this point if dmtcp checkpointing was NOT
    #requested (or this cluster doesn't support dmtcp checkpointing) and
    #the --phosim-checkpointing option was given when this program was run.
    #Defaults are no checkpointing of either type (in which case any
    #checkpointing , even if specified when running the phosim.py command is
    #ignored.
        
    #Determine if we are to enabble internal phosim checkpointing. This
    # can only happen if:
    # 1:DMTCP checkpointing was not requested on the command line
    # 2: AND phosim Internal checkpointing it was requested on the command line
    # 3: AND if phosim Internal checkpointing was requested on the 
    #    "./phosim -g condor" command line which was used to createthe dagman
    #    file.
    if  phosimChckPointReq:      
        #See if internal checkpointing was was instituted.
        internalCheckPointReq = isPhosimChckptRequested(inputDirFull)
    else:
        internalCheckPointReq = False
        phosimChckPointReq = False

    #Build the raytrace part of the file. For DMTCP or no checkpointing the
    #raytrace*_0 is all that we do. Only for the phosim internal
    #checkpointing do we invoke the rest of the raytrace*_n jobs. 
    if internalCheckPointReq:
         raytraceJobs = raytraceJobDict
    else:
         raytraceJobs =  raytraceJobDict_0
         

    # #################################################
    # At this point we make the raytrace jobs.
    # There are 4 possible scenarios here:
    # 1:--Trim-Disable==False:  If this script is run without the
    #    --Trim-Disable being given then the trim jobs are created and are
    #    submitted. The trimPIDs are placed in the jobNamePID dict. Then the
    #    raytrace jobs are created (replaceing any previously made raytrace
    #    jobs) using the trim dependencies.
    #    A: --Raytrace-Disable == False: If the --Raytrace-Disable option is
    #       also not given, the script will preceed to submit the newly
    #       created raytrace jobs.
    #    B: --Raytrace-Disable == True: If the  --Raytrace-disable is given
    #       the raytrace jobs will not be submitted (though they are created).
    # 2: --Trim-Disable == True: If this script is run  with the
    #    --Trim-Disable option given then the trim jobs will not be created
    #    or submitted. However any previously made raytrace jobs will still
    #    exist.
    #    A: --Raytrace-disable == False: If Raytrace-Disable is not given the
    #       existance of the raytrace jobs will be tested and if they exist
    #       they will be submitted. Otherwise the script will terminate with
    #       an error.
    #    B:  --Raytrace-disable == True: If Raytrace-Disable is also give
    #        nothing is done and an error is declared.


    # If the --E2adc-Disable option is chosen, then any time a raytrace job
    # is created the running of the e2adc task is not included.

    for jobName,raytraceSubmitFileName in raytraceJobs.items():

        #Get the important lines from the submission file
        jobSubDict = getSubmissionParams(raytraceSubmitFileName)

        #If internalCheckPointReq is False, then clear checkpointtotal and
        # checkpointcout to 0 every time by creating a new
        # .pars file and thus leave the original phosim -g condor produced
        # raytrace.pars file untouched in the work directory..
        # Otherwise leave the pars file alone.

        #Make a non internal checkpoint raytrace*.pars file if needed.
        parsFileName = jobSubDict['input']
        if not internalCheckPointReq:
            parsFileName = createNonCheckpointingParsFile(
                                                      jobSubDict['initialDir'], 
                                                      parsFileName)
            #print('parsFileName: ' + parsFileName)
            
        #Job cluster submission file. (Aribitrarily has .pbs extention even if
        #its a slrumn mediated cluster
        jobPBS = opt.workDir + "/" + jobName + '.pbs'

        # ##############################################################
        #At this point we need to see if the --Trim-Disable option was chosen
        #If it was and we get here, than we are assured the --Raytrace-Disable
        #option was not chosen and we are to try and submit a raytrace job.
        #If --Trim-Disable was not chosen we can just go ahead and make a new
        #raytrace job. If --Trim-Disable was chosen we do not create a new
        #raytrace job , but instead check to see if an old one exists and we
        #will submit it. If one does not exist, we declare a fatal error and
        #exit. 
        if not disableTrimSubmission:
            pfile=open(jobPBS,'w') #Don't forget to close this file when done
          
            #Build the cluster submission file.
            addClusterCommands(pfile, submitPBSList['INITIALLIST'],
                               submitPBSList['THREADCMD'], requestedThreads,
                               knl, haswellRegular)
            addDependancyCommand(pfile,jobName, dependancies, jobNamePID,
                                 submitPBSList['DEPENDCMD'])
            #At this point all batch commands have been added TOGETHER to pfile
            #now for cori, add the module load module liust commands
            if knl or haswellShared or haswellRegular :
                 addModuleLoadCommands(pfile,knl, haswellShared, haswellRegular)
                 
            addJobCommands(pfile, jobSubDict,'raytrace',jobName,
                           dmtcpChckPointReq, submitPBSList,parsFileName)

            # ##############################################################
            #Decide if this is the task we add the e2adc job to. 
            #1: If internal checkpointing we wait until we submit a task that
            #   is a parent to the e2adc child (last of checkpoint runs)
            #2: if no internal checkpointing than each chip only has 1
            #   raytrace*_0 task and we add e2adc to each jop.
            #get e2adcJobName
            e2adcJobName = 'e2adc_' + jobName[9:-2] 
            e2adcParentJobName =  dependancies[e2adcJobName]

            if internalCheckPointReq:
                e2adcParentJobName =  dependancies[e2adcJobName]
            else:
                #Puts e2adc into raytrace*_0 
                e2adcParentJobName = jobName 
             
            if  e2adcParentJobName == jobName:
                #Add the e2adc command. Get e2adc JobName from the raytrace 
                #JobName. 
                #Make sure we have this key in our e2adc dict and get the
                # e2adc submission file name (We could also check to see if
                # we had the dependency on the raytrace here in
                # dependencies also as a sanity check.)
                if (e2adcJobName not in e2adcJobDict) :
                    print('Key for e2adc job:' +  e2adcJobName + 
                          ' not in dagMan File' )
                    sys.exit(1)
                e2adcSubmitFileName = e2adcJobDict[e2adcJobName]
                e2adcSubDict = getSubmissionParams(e2adcSubmitFileName)
                e2adcParsFile = e2adcSubDict["input"]

                #This is where we we check to see if the 
                #--E2adc-Disable option was chosen
                if not disableE2adcSubmission: 
                    addTaskRunLine(pfile, e2adcSubDict, 'e2adc', e2adcParsFile,
                                   False, False)
                                             #will get executable path
                #addRemoveFileLine(pfile,e2adcSubDict['input'])
                                                           #Disabled for now

                #At this point we need to move the image fits files to the
                #output directory. First the completed "lsst_e" file name
                #from the .submit file transfer list. (other way to get
                #this name is to make it up from the jobName and to
                #search the e2adc*.pars file for the filter line!
                #This is easier.)

                transferList=e2adcSubDict['transferFiles'].split()
                for transFile in transferList:
                    if transFile[:6] == "lsst_e":
                        lsst_e_FileName = transFile
                        addMoveFileLine(pfile,lsst_e_FileName,
                                        opt.outputDir)
                        break
                        #Not really necessatry to get these _a_ files but
                        #its easy so do it. We can remove later(we did)
                        #if not useful. Use a wild card to get them all
                        #lsst_a_FileName = ("lsst_a" +
                        #                    lsst_e_FileName[6:-12] +
                        #                     '*' +".fits.gz")
                        #addMoveFileLine(pfile,lsst_a_FileName,
                        #                 opt.outputDir)
                        #break
                #The next line will probably work but lets wait a bit
                #before we enable it.
                #addRemoveFileLine(pfile,jobPBS) #Disabled for now

            #And thats all for the raytrace/e2adc .pbs file. Add date print and
            #close it!
            pfile.write('date' + '\n')
            pfile.close()

        #There are two conditions that must be met to submit this job.
        #1: The job exists
        #   and
        #2: The --Raytrace-Disable option was not chosen

        #now we can submit it. We may need the pid. 

        if not os.path.isfile(jobPBS):
             print("Fatal -- Attempt to submit raytrace job " + jobPBS )
             print('Fatal -- Job does not exists. Check "Disable" options!')
             sys.exit(1)
             
        if not disableRaytraceSubmission:
 
            submitCmd = ( submitPBSList['SUBMITCMD'] + 
                          ' -e ' + jobSubDict['initialDir'] + '/errors/' + 
                          jobName + '.pbs.err ' + 
                          ' -o ' + jobSubDict['initialDir']  +  '/logs/' + 
                          jobName + 'pbs.log ' +   jobPBS)
  
            #print ('Submit cmd: ' + submitCmd) 

            # #################################################
            subPBSLine = subprocess.check_output(submitCmd, shell=True)
        
            #Debug lines follows
            #jobCount = jobCount +1
            #subPBSLine = generateDebugSubPBSLine(jobCount)

            print(jobName + ': ' + subPBSLine)


            #Pick up the PID for dependencies later
            # We should probably check the return code here. Later!
            JobPid = getPid(jobPBS,subPBSLine)
            # #################################################
            #save this in a dict with the job name as a key. Used
            # when setting up internal checkpoint raytrace dependencies
            jobNamePID[jobName]=JobPid

            #We should probably check the status here (did command succdeed?)
            #Later!

            #cleanup
	        #removeFile(raytraceSubmitFileName)
            #removeFile(e2adcSubmitFileName)

    return

##main function
# First we will only be interested in getting the trim jobs running. We will
# save the trim  submission job ID's for later use as a dependency requirment
# for the raytrace jobs
def main():
     print('cluster_submit.py: v1.3')

     # Get the current directories as defaults. This sssumes we are running
     # from the tools direcoty. These may all be overridden by commandl line
     # arguments
     toolsDir=os.path.split(os.path.abspath(__file__))[0]
     phosimDir = (toolsDir + "/..")
     defaultOutputDir=os.getenv('PHOSIM_OUTPUT_DIR', phosimDir+'/output')
     defaultWorkDir=os.getenv('PHOSIM_WORK_DIR', phosimDir+'/work')

     # Parse the command line
     parser = optparse.OptionParser(usage='%prog dagManFile ' +
	                            ' [<arg1> <arg2> ...]')

     # define acceptable options
     parser.add_option('-w','--work',dest="workDir",default=defaultWorkDir,
                       help='temporary work directory')
     parser.add_option('-o','--output',dest="outputDir",
                       default=defaultOutputDir, help='output directory')
     parser.add_option("--Trim-Disable", action="store_true",
                       dest="disableTrim", default=False,
                       help= ("Disable creation and submission of the new "
                              "Trim batch jobs. New raytrace/e2adc jobs are "
                              "not created but previously made raytrace/e2adc "
                              "jobs with previous trim job dependancies may "
                              "be submitted (but see --Raytrace-disable). "
                              "Default [%default] is to submit trim jobs and "
                              "create new raytrace/e2adc jobs (with the new "
                              "trim jobs as dependancies).") )
     parser.add_option("--Raytrace-Disable", action="store_true",
                       dest="disableRaytrace", default=False,
                       help= ("Disable submission of the Raytrace/e2adc batch "
                              "jobs. Default [%default] is to submit the "
                              "raytrace/e2adc jobs (if they exist, see "
                              "--Trim-Disable).") )
     parser.add_option("--E2adc-Disable", action="store_true",
                       dest="disableE2adc", default=False,
                       help= ("Disable inclusion of running of the e2adc task "
                              "within any generated raytrace job. Default " +
                              "[default: %default] is to include the running "
                              "of e2adc tasks in raytrace jobs." ) )
     parser.add_option("--dmtcp-checkpoint", action="store_true",
                       dest="dmtcpChckpntReq", default=False,
                       help= ( "setup checkpoint of raytrace tasks using " +
                               "dmtcp. Internal phosim checkpointng ignored" +
                               "[default: %default]") )
     parser.add_option("--phosim-checkpoint", action="store_true",
                       dest="phosimChckpntReq", default=False,
                       help= ("Use internal phosim checkpointing scheme if " +
                              "checkpointtotal > 0 in raytrace*.pars files " +
                              "[default: %default]" + "(Mainly for use on "
                              "Purdue clusters, use dmtcp checkpointing on "
                              "others.)" ) )
     parser.add_option("--Cori-KNL", action="store_true",
                       dest="coriKNLPartition", default=False,
                       help= ("If running on the Cori cluster, define partition " +
                              "to use to be the KnightsLanding partition.)" ) )
     parser.add_option("--Cori-Haswell-Shared", action="store_true",
                       dest="coriHaswellSharedPartition", default=False,
                       help= ("If running on the Cori cluster, define partition " +
                              "to use to be the Haswell partition and Shared "
                              "queue (default for Cori)" ) )
     parser.add_option("--Cori-Haswell-Regular", action="store_true",
                       dest="coriHaswellRegularPartition", default=False,
                       help= ("If running on the Cori cluster, define partition " +
                              "to use to be the Haswell partition and Regular "
                              "queue (default for Cori)" ) )
     parser.add_option("--Cori-Phosim-Version", 
                       dest="PhosimVersion", default=phosimVersion,
                       help= ("If running on the Cori cluster, define phosim "
                              "version to use to be used (EX: dev, 3.7beta, 3.6.1 "
                              "default=dev)" ) )
     if len(sys.argv)<2:
          usage()
          sys.exit()
     if sys.argv[1] in ('-h', '--help'):
          usage()
          sys.exit()
               
     # parse_args returns a pair of values, remainder should have dagMan
     # file name
     opt, remainder = parser.parse_args(sys.argv[1:]) 

     #require dagMan file path/name
     if len(remainder) == 0:
          print 'dagMan File name argument is required'
          usage()
          sys.exit()
     if len(remainder) != 1:
          print ('Too many arguments on command line: %s',  remainder)
          usage()
          sys.exit()
     
     dagManFile = remainder[0]
     # see if we have a path to this file or if we should assume the work
     #directory
     dagDir=os.path.split(dagManFile)[0]
     if len(dagDir) == 0:
          dagManFileFull = opt.workDir + '/' + dagManFile
     else:
          dagManFileFull = dagManFile
     
     #Probably should Make sure paths to work/logs and work/errors

     
     createAndSubmitJobs(opt, dagManFileFull)
          
if __name__ == "__main__":
    main()
