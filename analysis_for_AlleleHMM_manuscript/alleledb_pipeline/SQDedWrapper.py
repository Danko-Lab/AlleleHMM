#!/usr/bin/env /home2/njc2/myInstalls/bin/njcPython
# Set the above path as part of the install?

import logging, nws.client, os, signal, socket, sys, thread, time
LogC, LogD, LogE, LogI, LogW = logging.critical, logging.debug, logging.error, logging.info, logging.warning

myHost = socket.gethostname()

prog, verbose, logFilePath, key, nwssHost, nwssPort, jobIndex = sys.argv

logging.basicConfig(filename='%s%swrapper_%s.log'%(logFilePath, os.sep, jobIndex), format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
LogI('Wrapper script running on %s for job %s (pid: %d).'%(myHost, jobIndex, os.getpid()))

verbose = verbose == 'True'

# Used by the sentinel thread to check if the runIt thread's wait has returned.
childExitedLock = thread.allocate_lock()
childExitedLock.acquire()

# Used by the sentinel thread to wait for the child to be forked.
childExistsLock = thread.allocate_lock()
childExistsLock.acquire()

childPid = None

# Used by the main thread to wait for the sentinel thread to exit.
sentinelDoneLock = thread.allocate_lock()
sentinelDoneLock.acquire()

# Used by the main thread to wait for the run thread to exit.
runItDoneLock = thread.allocate_lock()
runItDoneLock.acquire()

KILLWAITDECIS = 30 # Wait (in deciseconds) for KILL SIGNAL to have an effect.
TERMWAITDECIS = 60 # Wait (in deciseconds) for TERM SIGNAL to have an effect.

testing = False
if testing:
    # Referencing the workspace here ensures the main thread's connection owns the workspace.
    ws = nws.client.NetWorkSpace(key, serverHost=nwssHost, serverPort=int(nwssPort))

def sentinel():
    # Make sure that the child was instantiated.
    childExistsLock.acquire()

    ws = nws.client.NetWorkSpace(key, serverHost=nwssHost, serverPort=int(nwssPort))
    sdc = -1
    try: sdc = ws.find('shut down %s'%jobIndex)
    except Exception, e:
        if not childExitedLock.locked():
            LogI('Sentinel failed in find, but child has exited.')
            sentinelDoneLock.release()
            return
        LogI('Sentinel failed in find, attempting shut down.')
        
    LogI('Sentinel got shut down code %d.'%sdc)
    if sdc == 0:
        # We sent this to ourselves, go quietly.
        sentinelDoneLock.release()
        return

    LogI('Sentinel sending SIGTERM to process group %d.'%childPid)
    try: os.killpg(childPid, signal.SIGTERM)
    except OSError, e:
        if e.errno == 3: # The child is gone.
            sentinelDoneLock.release()
            return
        LogI('Problem with sending SIGTERM: %s.'%str(e))

    for x in range(TERMWAITDECIS):
        if not childExitedLock.locked(): break
        time.sleep(.1)
    else:
        LogI('Sentinel sending SIGKILL to process group %d.'%childPid)
        try: os.killpg(childPid, signal.SIGKILL)
        except OSError, e:
            if e.errno == 3: # The child is gone.
                sentinelDoneLock.release()
                return
        LogI('Problem with sending SIGKILL: %s.'%str(e))

        for x in range(KILLWAITDECIS):
            if not childExitedLock.locked(): break
            time.sleep(.1)
        else:
            if childExitedLock.locked():
                LogI('Looks like %d won\'t die.'%childPid)
                ws.store('Job %s Status'%jobIndex, (None, startTime, time.time(), 1, childPid, myHost))
                runItDoneLock.release() # Do this on behalf of the runIt thread, which will be stuck waiting for the child.
    sentinelDoneLock.release()
thread.start_new_thread(sentinel, ())

def runIt():
    global childPid

    ws = nws.client.NetWorkSpace(key, serverHost=nwssHost, serverPort=int(nwssPort))
    LogI('Opened workspace %s on %s at %s.'%(key, nwssHost, nwssPort))

    jobTag = 'job %s'%jobIndex
    cmd = ws.fetch(jobTag)
    LogI('Fetched job: %s.'%cmd)
    try:    ws.deleteVar(jobTag)
    except: pass

    startTime = time.time()

    ecmd = ['/bin/bash', '-c', cmd]
    childPid = os.fork()
    if childPid == 0:
        # Create a new process group to simplify killing the child and its descendants.
        os.setpgrp()
        os.execv(ecmd[0], ecmd)

    LogI('Child %d started on %s at %f.'%(childPid, myHost, startTime))
    childExistsLock.release()

    while 1:
        (retpid, retval) = os.waitpid(childPid, 0)
        if retpid: break

    LogI('Child %d returned on %d.'%(childPid, retval))
    childExitedLock.release()
    ws.store('Job %s Status'%jobIndex, (retval, startTime, time.time(), 0, childPid, myHost))
    ws.store('shut down %s'%jobIndex, 0) # let the sentinel know that the child is done.
        
    runItDoneLock.release()
thread.start_new_thread(runIt, ())

if testing:
    LogI('Storing job tuple.')
    ws.store('job %s'%jobIndex, './testJob')
    LogI('Sleeping.')
    time.sleep(7)
    LogI('Storing shut down tuple.')
    ws.store('shut down %s'%jobIndex, 1)

# Wait for the threads to finish.
runItDoneLock.acquire()
sentinelDoneLock.acquire()

if testing:
    print ws.fetch('Job %s Status'%jobIndex)

LogI('Wrapper script finished on %s for job %s (pid: %d).'%(myHost, jobIndex, os.getpid()))
sys.exit(0)
