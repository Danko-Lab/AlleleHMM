#!/usr/bin/env /home2/njc2/myInstalls/bin/njcPython

# wrapper log to run/job specific file (based on key?)
# => add log a la sqDed to wrapper

import logging, nws, os, pnwsst, Queue, signal, sys, thread, threading, time
LogC, LogD, LogE, LogI, LogW = logging.critical, logging.debug, logging.error, logging.info, logging.warning

# we assume that the helper script lives in the same directory as this script.
WRAPPER = os.path.dirname(os.path.realpath(__file__))+os.path.sep+'SQDedWrapper.py'

#import simpleQueueConfig

# This variation assumes that the nodes are dedicated to this run
# (e.g., via a PBS allocation), so we don't need to check loads on
# the nodes. We just need to track assignments and their completions.

# A set of nodefiles is provided.  By default, the program will look for
# $PBS_NODEFILE.  This can be prevented via the --nopbs flag, and
# additional nodefiles can be provided on the command line.

# These files will be watched, and if their timestamps change, the list
# of nodes to be used will be reloaded.  Running jobs will not be
# stopped, however.

# If the job is terminated prematurely, or if any jobs fail, a new joblist file
# will be created (using the original filename with '.REMAINING' appended) that
# contains all jobs that didn't apparently finish successfully.

class JobList:
    done, jobs = {}, []

    def __init__(self, f):
        self.jobFileName = f
        for l in open(f):
            l = l.strip()
            if not l or l[0] == '#': continue
            self.jobs.append(l)

    def setDone(self, i): self.done[i] = True

    def dumpRemaining(self):
        # doing this even when all jobs are "done" avoids stale REMAINING files.
        dumpfile=self.jobFileName+'.REMAINING'
        f = open(dumpfile, 'w')
        if len(self.jobs) > len(self.done):
            LogI('Remaining jobs dumped to %s.'%dumpfile)
            for i in xrange(len(self.jobs)):
                if i not in self.done:
                    f.write("%s\n" % self.jobs[i])
        f.close()
        
pbsFile = None
class NodePool:
    def __init__(self, key, nodeFiles, ws):
        self.activeJobs = {}
        self.fileList = []
        self.key = key
        self.ws = ws
        
        LogI('Processing node files: '+str(nodeFiles))
        for f in nodeFiles:
            if f == '$PBS_NODEFILE':
                if 'PBS_NODEFILE' in os.environ:
                    global pbsFile
                    pbsFile=os.environ['PBS_NODEFILE']
                    self.fileList.append([pbsFile, os.stat(pbsFile)[8]])
                else:
                    LogE('Error, PBS_NODEFILE not defined')
                    sys.exit(1)
            else:
                self.fileList.append([f, os.stat(f)[8]])

        LogD('Filelist: %s'%self.fileList)

        self.drainedSem = threading.Semaphore(0)
        self.draining = False

        self.nodeQueue = Queue.Queue(0)
        self.nodes = {}
        self.monSem = threading.Semaphore()       # used to implement monitor-like access to an instance.
        self.resetNodePool()

    def resetNodePool(self):
        # monitor is acquired by checkNodeFiles
        # first set all allowed jobs to 0
        for r in self.nodes.itervalues():
            r[0] = 0
        # go through all node files, and increment allowed jobs for each entry
        for f, time in self.fileList:
            for n in open(f):
                n = n.strip()
                if not n or n[0]=='#': continue
                r = self.nodes.get(n, [0, 0, []]) # per node data: allowed jobs, assigned jobs, extra info (unused)
                r[0] += 1
                self.nodes[n] = r

        # apply PBS limit. doing it here takes a conservative approach
        # to the case where a node is listed in both the PBS nodefile
        # and another file.
        if pbsFile:
            for n in open(pbsFile):
                n = n.strip()
                if not n or n[0]=='#': continue
                self.nodes[n][0] = min(self.nodes[n][0], oArgs.maxJobsPerNode)

        # rebuild queue
        while not self.nodeQueue.empty(): self.nodeQueue.get()
        for n, r in self.nodes.iteritems():
            for x in xrange(max(0, r[0]-r[1])):
                LogI('Adding %s to node queue.'%n)
                self.nodeQueue.put(n)

        LogI('Nodelist: %s' % self.nodes.keys())

    def checkNodeFiles(self):
        self.monSem.acquire()
        changed=False
        for e in self.fileList:
            current=os.stat(e[0])[8]
            if current > e[1]:
                LogI('Nodefile %s changed'%e[0])
                e[1]=current
                changed=True
        if changed: self.resetNodePool()
        self.monSem.release()
            
    def drainPool(self):
        # wait for all assigned jobs to finish.
        self.monSem.acquire()
        outstanding = sum([r[1] for r in self.nodes.itervalues()])
        self.draining = True
        self.monSem.release()                
        if outstanding: self.drainedSem.acquire()

    def getNode(self, jobIndex):
        # return a node that can be assigned a new job --- block if
        # necessary for one to become available.
        n = self.nodeQueue.get(True)
        if n == 'shutting down':
            self.nodeQueue.put(n)
        else:
            self.monSem.acquire()
            self.activeJobs[jobIndex] = True
            self.nodes[n][1] += 1
            self.monSem.release()

        return n

    def releaseNode(self, n, jl, i, completed=True):
        LogI('Releasing %s for job %d (%s).'%(n, i, completed))
        self.monSem.acquire()
        if self.activeJobs: self.activeJobs.pop(i)
        if completed: jl.setDone(i)
        r = self.nodes[n]
        r[1] -= 1
        self.nodes[n] = r
        if self.draining:
            if not sum([r[1] for r in self.nodes.itervalues()]):
                LogI('Releasing drainedSem.')
                self.drainedSem.release()
        else:
            if r[1] < r[0]: self.nodeQueue.put(n)
        self.monSem.release()

    def shutdown(self):
        # this should only be called by the signal handler.
        LogI('In shutdown.')
        self.monSem.acquire()
        # empty node queue.
        while not self.nodeQueue.empty(): self.nodeQueue.get()
        LogI('Node queue is now empty.')
        self.nodeQueue.put('shutting down')
        self.draining = True
        self.jobsToShutdown = self.activeJobs
        self.activeJobs = None # Shouldn't try to use this as a dictionary from here on out.
        self.monSem.release()
        LogI('Poison node in queue.')
        
def runJob(np, jl, n, i, j, key):
    wsSem.acquire()
    ws.store('Launched', str(int(ws.fetch('Launched'))+1))
    wsSem.release()
    logFilePath = os.path.dirname(os.path.realpath(oArgs.logFile)) # all other log files are placed in the same directory as the main log file.
    if oArgs.noWrapper:
        cmd = 'ssh -n %s %s'%(n, j)
    else:
        j = j.strip()
        if j[0] in '\'"':
            LogW('removing quotes from: '+j)
            j = j[1:-1]
        wsSem.acquire()
        ws.store('job %d'%i, j)
        wsSem.release()
        cmd = 'ssh -nf %s "%s %s %s %s %s %d %d 2> %s/%d.ssh.err"'%(n, WRAPPER, oArgs.wrapperVerbose, logFilePath, key, oArgs.nwssHost, oArgs.nwssPort, i, logFilePath, i)
    # is logging thread safe? yes.
    LogI('Launching:\n\t%s\n\t%s' % (cmd, j))
    status = os.system(cmd)
    LogI('Job %d launch reported %d' % (i, status))
        
    # be careful in the following not to block or blow up.
    wsSem.acquire()
    try:
        me = 'Job %d Status'%i
        LogI('Fetching extra data for job %d from workspace.'%i)
        while 1:
            status, startTime, stopTime, rogue, pid, host = ws.fetchTry(me, (-1, -1., -1., True, -1, 'no where'))
            if pid != -1: break
            wsSem.release()
            time.sleep(1)
            wsSem.acquire()
        try:    ws.deleteVar(me)
        except: pass
        if status != None: status = divmod(status, 256)
        LogI('Extra data for job %d from workspace: %s %f %f %s %d %s'%(i, status, startTime, stopTime, rogue, pid, host))
        jStatus.write('%d\t%s\t%f\t%f\t%d\t%d\t%s\t%s\n'%(i, status, startTime, stopTime, rogue, pid, host, j))
        if rogue:
            jRogues.write('%s\t%d\n'%(host, pid))
        ws.store('Done', str(int(ws.fetch('Done')) + 1))
        if status == 0:
            ws.store('Succeeded', str(int(ws.fetch('Succeeded')) + 1))
        else:
            ws.store('Failed', str(int(ws.fetch('Failed')) + 1))
    except Exception, e:
        LogI('Encountered exception "%s" while accessing extra status info for job %d.'%(e, i))
    wsSem.release()
        
    if status == (0, 0):
        # Everything OK.
        np.releaseNode(n, jl, i)
    elif status and status[0] == 0 and oArgs.ignoreErrors:
        # Means the job exited with a non-zero exit code, but we've been told to ignore that.
        np.releaseNode(n, jl, i)
    else:
        # Either the job was terminated, it exited with a non-zero
        # code which we are not ignoring, or it may be a rogue. Flag
        # as incomplete.
        np.releaseNode(n, jl, i, completed=False)

def runJobs(np, jl, key):
    numJobs = len(jl.jobs)
    ws.store('Total jobs', str(numJobs))
    LogI('Total jobs %d.'%numJobs)

    allLaunched = False
    for i in xrange(len(jl.jobs)):
        # must call getnode in this context *before* thread is created.
        n = np.getNode(i)
        if n == 'shutting down': break
        thread.start_new_thread(runJob, (np, jl, n, i, jl.jobs[i], key))
    else:
        LogI('All jobs launched.')
        allLaunched = True

    if not allLaunched: LogI('Looks like we\'re shutting down, so skipping jobs from index %d on.'%i)
    LogI('Draining pool.')
    np.drainPool()

def setupLogging(logLevel, logFile):
    logging.basicConfig(level=logLevel,
                        filename=logFile,
                        filemode='w',
                        format='%(asctime)s.%(msecs)03d %(levelname)s %(message)s',
                        datefmt='%Y-%m-%d %H:%M:%S')
                                           
    console = logging.StreamHandler()
    console.setLevel(logging.INFO)
    formatter=logging.Formatter('%(levelname)s %(asctime)s %(message)s')
    console.setFormatter(formatter)
    logging.getLogger('').addHandler(console)


if '__main__' == __name__:
    import optparse

    nwssSet = False
    def monitor_store(option, opt_str, value, parser):
        global nwssSet
        nwssSet = True
        setattr(parser.values, option.dest, value)

    opts = optparse.OptionParser(description='Process a simple queue of jobs using a dedicated collection of nodes.',
                                 usage='usage: %prog [options] JobFile')
    opts.add_option('-H', '--nwssHost', help='nws server host (%default).', action='callback', callback=monitor_store, metavar='Host', default='localhost') 
    opts.add_option('-i', '--ignoreErrors', help='Consider a job done even if it returns an error code.', action='store_true', default=False)
    opts.add_option('-l', '--logFile', help='Specify log file name (%default).', metavar='LogFile', default='log.out') 
    opts.add_option('--maxJobsPerNode', help='When running with PBS, limit jobs to N per node.',  metavar='N', type='int', default=1000000) 
    opts.add_option('-n', '--nodeFiles', help='Comma separated list of files listing nodes to use (%default).', metavar='FileList', default='$PBS_NODEFILE')
    opts.add_option('--noWrapper', help='Do not use a wrapper process to monitor task execution.', action='store_true', default=False)
    opts.add_option('-p', '--nwssPort', help='nws server port (%default).', action='callback', callback=monitor_store, metavar='Port', type='int', default='8765') 
    opts.add_option('--pnwss', help='Create a personal workspace for this run.', action='store_true', default=False)
    opts.add_option('-v', '--verbose', action='store_const', const=logging.DEBUG, default=logging.INFO)
    opts.add_option('-V', '--wrapperVerbose', action='store_true', default=False)
       
    oArgs, pArgs = opts.parse_args()
    if len(pArgs) != 1:
        # Eventually change to allow multiple job files, with some sort of indication of which file a job comes from.
        print >>sys.stderr, 'Need one (and only one) job file.'
        opts.print_help(sys.stderr)
        sys.exit(1)

    setupLogging(oArgs.verbose, oArgs.logFile)
    if oArgs.pnwss:
        if nwssSet: LogW('pnwss overrides nwssHost and nwssPort.')
        pnwss = pnwsst.NwsLocalServer(logFile=oArgs.logFile+'.pnwss')
        LogI('Started personal networkspace server: %s %d %d'%(pnwss.host, pnwss.port, pnwss.webport))
        setattr(oArgs, 'nwssHost', pnwss.host)
        setattr(oArgs, 'nwssPort', pnwss.port)

    LogI('Control process is %d.'%os.getpid())
    LogI('sqDedicated run started using nws server %s %d'%(oArgs.nwssHost, oArgs.nwssPort))
    
    key=('SENTINEL %s %s %s' % (os.environ['LOGNAME'], time.asctime(), pArgs[0])).replace(' ', '_')
    LogI('Sentinel key: %s' % key)

    ws = nws.client.NetWorkSpace(key, serverHost=oArgs.nwssHost, serverPort=oArgs.nwssPort)
    ws.store('Launched', '0')
    ws.store('Done', '0')
    ws.store('Failed', '0')
    ws.store('Succeeded', '0')

    # needed to synchronize access by threads running jobs.
    wsSem = threading.Semaphore()

    np = NodePool(key, oArgs.nodeFiles.split(','), ws)
    jl = JobList(pArgs[0])

    jStatus = open(jl.jobFileName+'.STATUS', 'w', 0)
    jRogues = open(jl.jobFileName+'.ROGUES', 'w', 0)

    goodbye = False

    def runThread():
        runJobs(np, jl, key)
        LogI('runJobs has finished.')
        jRogues.close()
        jStatus.close()
        
        jl.dumpRemaining()
        LogI('Run completed.')
        pnwss.stop()
        global goodbye
        goodbye = True
        os.kill(os.getpid(), signal.SIGTERM)

    rt = threading.Thread(None, runThread)
    rt.start()

    def noteSignal(sig, frame):
        LogI('Ignoring signal %d.'%sig)
        
    for x in range(1, signal.NSIG):
        if x == signal.SIGCHLD: continue
        try: signal.signal(x, noteSignal)
        except: pass

    def shutItDown(sig, frame):
        # add signal.alarm to as a deadman switch?
        if goodbye:
            LogI('Received signal %d, but already shutting down.'%sig)
            return
        else:
            LogI('Received signal %d. Shutting down.'%sig)
        np.shutdown()
        wsSem.acquire()
        for i in np.jobsToShutdown:
            ws.store('shut down %d'%i, 1)
            LogI('Told job %d to shutdown.'%i)
        wsSem.release()
        
    sigs = [ 1, 2, 3, 15 ]
    for x in sigs:
        try: signal.signal(x, shutItDown)
        except: pass

    while 1:
        signal.pause()
        if goodbye: break

    LogI('Falling off the end of the world.')
