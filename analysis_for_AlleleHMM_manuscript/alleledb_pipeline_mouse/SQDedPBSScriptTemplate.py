#!/usr/bin/env /home2/njc2/myInstalls/bin/njcPython
#PBS -q %(queue)s %(group)s -A %(account)s -N %(title)s
#PBS -l nodes=%(numProcs)s:ppn=4
#PBS -o PBS_%(title)s_out.txt
#PBS -e PBS_%(title)s_err.txt

import logging, os, signal, sys, time
LogC, LogD, LogE, LogI, LogW = logging.critical, logging.debug, logging.error, logging.info, logging.warning

os.chdir(os.getenv('PBS_O_WORKDIR'))
sqDir = 'SQ_Files_'+os.getenv('PBS_JOBID')+os.sep
os.mkdir(sqDir)

# Set a few initial arguments for a function.
def setArgs(func, *args):
    def invokeFunc(*extras):
        return func(*(args+extras))
    return invokeFunc

# Establish func as the handler for each signal in sigs. 
def setupHandlers(sigs, func):
    for x in sigs:
        try:    signal.signal(x, func)
        except RuntimeError, e:
            if e.args[0] == 22:
                pass # Ignore this---we tried to set the handler for a signal that cannot be caught.
            else:
                LogE('While setting a signal handler for %%d: %%s.'%%(x, str(e)))
            
# A generic handler to relay a signal to a pid.
def relaySignal(pid, logIt, sig, frame):
    if logIt: LogI('Sending %%d to %%d.'%%(sig, pid))
    os.kill(pid, sig)

# pbs_mom aggressively kills processes under its watch. Do an initial
# fork to create a child that will be in a different session and so
# will be able to shut down in a more orderly fashion.
firstPid = os.fork()
if firstPid:
    logging.basicConfig(filename=sqDir+'PBS_script_0.log', format='%%(asctime)s %%(levelname)s %%(message)s', level=logging.DEBUG)
    LogI('Batch script starting (pid: %%d).'%%os.getpid())
    # This is the parent process.  Set up a handler for every signal
    # except SIGCHLD. The handler relays the signal to the child
    # process identified by firstPid.
    setupHandlers([x for x in range(1, signal.NSIG) if x != signal.SIGCHLD],
                  setArgs(relaySignal, firstPid, False))
    r = (-1, -1)
    while 1:
        # do nothing but wait for the newly created child process to exit.
        try: 
            r = os.waitpid(firstPid, 0)
            if r[0]: break
        except OSError, e:
            if e.errno != 4:
                LogE(str(e))
                break
    LogI('Batch script ending (%%s).'%%str(r))
    sys.exit(0)

# Create a new session to elude the gaze of pbs_mom's watchful eye.
os.setsid()

# Now we can set up and run the simple queue driver script.

logging.basicConfig(filename=sqDir+'PBS_script_1.log', format='%%(asctime)s %%(levelname)s %%(message)s', level=logging.DEBUG)
LogI('Batch script starting in earnest (pid: %%d).'%%os.getpid())

# Command to invoke the simple queue driver script
cmd = ['%(sqScript)s',
       '--logFile=%%s/SQ.log'%%sqDir,
       '--maxJobsPerNode=%(mppn)s',
       '--pnwss',
       '--wrapperVerbose',
       '%(jobFile)s']
LogI('About to fork/exec %%s.'%%cmd)

# Fork/exec a child to run the script.
secondPid =  os.fork()
if secondPid == 0: os.execv(cmd[0], cmd)
    
LogI('%%d forked %%d.'%%(os.getpid(), secondPid))

setupHandlers([x for x in range(1, signal.NSIG) if x != signal.SIGCHLD],
              setArgs(relaySignal, secondPid, True))

# Now wait for the simple queue driver to exit.
r = (-1, -1)
while 1:
    try: 
        LogI('Waiting for %%d.'%%secondPid)
        r = os.waitpid(secondPid, 0)
        LogI('waitpid returned: %%s.'%%str(r))
        if r[0]: break
    except OSError, e:
        if e.errno != 4:
            LogE(str(e))
            break

LogI('Batch script exiting, %%s.'%%str(r))
