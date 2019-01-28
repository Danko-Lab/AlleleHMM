#!/usr/bin/env /home2/njc2/myInstalls/bin/njcPython

# args q, numProcs[.JobsPerProc], account, title, job file
import optparse, os, sys

qs = {
'amm': ('molinaro', 'molinaro'),
'eph': ('eph', 'eph'),
'fs': ('sigworth', 'sigworth'),
'gen': ('general', ''),
'mg': ('gerstein', 'gerstein'),
'd': ('default', ""),
}
nns = qs.keys()
nns.sort()
nns = ', '.join(nns)

opts = optparse.OptionParser(usage='''%%prog queueInfo NumNode[.MaxProcsPerNode] Account Title TaskFile

Generate a PBS submission script to run a job that works through the
simple queue of tasks given in the TaskFile (each task is a quoted
command line). The job will distribute the tasks one by one to the
processors allocated by PBS, sending a new task to a processor when it
has completed a previous one. The jobs exits when all tasks have
completed. Various logging files contain information about the
execution of individual tasks and the overall state of the task
processing. Files named TaskFile.XXXX report summary data.

The file pnwss_XXXX.info provides contact information (host: first
field, port: third field) for a web interface to track the overall
status (total tasks, done, succeeded, failed).


NOTE: You must submit the generated script to PBS to actually run the job.

queueInfo may be one of the abbreviations:

   %s

or queueName[:group]. If the latter and group is not specified it is
assumed to be the same as the queueName.

MaxProcsPerNode defaults to 4 if not given.'''%nns)

oArgs, pArgs = opts.parse_args()
if len(pArgs) != 5:
    opts.print_help(sys.stderr)
    sys.exit(1)

q, numProcs, account, title, jobFile = pArgs

if q in qs:
    queue, group = qs[q]
else:
    if ':' in q:
        queue, group = q.split(':', 1)
    else:
        queue, group = q, q

if queue == 'general':
    print >>sys.stderr('''\
Keep in mind when submitting to queue "general" that jobs running on
this queue can be pre-empted.  An effort is made to clean up, but it
would be wise to check for run-away processes if this submission does
not run to completion.''')

if group: group = '-W group_list=%s'%group

mppn = '4'
if '.' in numProcs:
    numProcs, mppn = numProcs.split('.', 1)

# We assume that related script lives in the same directory as this script.
myDir = os.path.dirname(os.path.realpath(__file__))+os.path.sep
sqScript = myDir + 'SQDedDriver.py'
PBSScript = open(myDir + 'SQDedPBSScriptTemplate.py').read()%locals()

print PBSScript
