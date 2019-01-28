
'''
This script blocks until the file(s) passed as arguments all exist.  It's useful for SimpleQueue
'''

import sys, time, os.path

waitfor = sys.argv[1:]

while True:
    wait=False
    for f in waitfor:
        if not os.path.exists(f): 
            wait=True
            break
    if wait==False:
        sys.exit(0)
    time.sleep(10)
