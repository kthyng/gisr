"""
Attempt at having a control-all file for oil modeling paper.

Run commands are included in this file, so to run, just do
`python run.py` from the command prompt and include flags with
-- in front to run whatever subprojects you want. Options include:
sensitivity, outer_f, galv_b, bara_b, and dwh_f for projects, and
compile to run the latex compilation.
So, to run several subprojects and compile, type
`python run.py --sensitivity --galv_b --compile` in the terminal.
Also note that this is meant to run in python version 2.7.
"""

import matplotlib as mpl
mpl.use("Agg") # set matplotlib to use the backend that does not require a windowing system
import numpy as np
import sys
import os
import netCDF4 as netCDF
import pdb
import matplotlib.pyplot as plt
import glob
from datetime import datetime, timedelta
import time
import tracpy
import init
import projects
import subprocess

# Units for time conversion with netCDF.num2date and .date2num
units = 'seconds since 1970-01-01'

if __name__ == "__main__":

    # Make sure necessary directories exist
    if not os.path.exists('logs'):
        os.makedirs('logs')

    process_queue = []
    poll_interval = 600.0 # number of seconds to wait

    # List of possible subprojects to run
    arg_switches = {"sensitivity":False, "outer_f":False, "galv_b":False,
                    "bara_b":False, "dwh_f":False, "galv_f":False}

    # Check to see if any input arguments were on the command line
    if len(sys.argv) > 1:
        for arg in sys.argv[1:]:
            # Select out all of the flags sent in
            # Flags should be indicated using '--' in front
            key = arg[2:] # take off preceding '--' from flags

            # Set all input keys to true
            if key in arg_switches.keys():
                arg_switches[key] = True
            else:
                raise ValueError("Unknown test %s requested." % key)

    # Find number of cores to use
    # Use the environmental variable if set
    # this is how one can override the computer setting if desired
    if os.environ.has_key('OMP_NUM_THREADS'):
        max_processes = int(os.environ['OMP_NUM_THREADS'])
    # Otherwise find out the number of cores on the computer
    else:
        if 'Linux' in os.uname(): # Linux
            # Find number of cores on a Linux box
            temp = os.popen('less /proc/cpuinfo | grep processor').read()
            max_processes = temp.count('processor')

        elif 'Darwin' in os.uname(): # Mac
            # Find number of cores on a mac
            max_processes = os.popen('system_profiler | grep "Cores"').read()[-2]

    # Make list of how many processes can run
    proc_count = list(np.arange(1,max_processes+1))

    # Prepare date for log file
    tm = time.localtime()
    year = str(tm[0]).zfill(4)
    month = str(tm[1]).zfill(2)
    day = str(tm[2]).zfill(2)
    hour = str(tm[3]).zfill(2)
    minute = str(tm[4]).zfill(2)
    second = str(tm[5]).zfill(2)
    date = '%s-%s-%s-%s:%s.%s' % (year,month,day,hour,minute,second)

    cmd_list = []
    log_name = []
    # Compile calls to do into a list so they can be distributed
    # amongst available cores 
    # On a Linux machine, the call assigns the process to a specific 
    # core by number using `taskset [number]`. Can't do this on a mac,
    # but hopefully the computer will distribute properly.
    # On all machine, the commands will redirect the screen output
    # to a log file, and the process will run in the background
    for (key,value) in arg_switches.iteritems():
        if value:
            # Set up log file
            log_name.append('%s-%s' % (date,key))

            # # Set up command
            # cmd_list.append('python2.7 projects/%s.py' % key)
            # Set up command
            cmd_list.append('python2.7 %s.py' % key)
                            
            # if this is a linux machine, keep track of which 
            # core we are using
            if 'Linux' in os.uname(): # Linux
                # Control the core to run on if on a Linux machine
                cmd_list[-1] = 'taskset %s %s' % (proc_count[0],cmd_list[-1])
                # move first core number to end of list
                proc_count.append(proc_count.pop(0))

            # when the max # of processes are compiled and running
            # check periodically to see if the next should be added
            while len(process_queue) == max_processes:
                print "Number of processes currently:",len(process_queue)
                for process in process_queue:
                    if process.poll() == 0:
                        # print 'process ended:' process
                        process_queue.remove(process)
                time.sleep(poll_interval)

            # Add on process to queue list
            # Remove the used log file from the list
            log_file = open('logs/' + log_name.pop(0) + '.txt','w')
            log_file.write('Started ' + date)
            print cmd_list[0]
            print log_file
            process_queue.append(subprocess.Popen(cmd_list.pop(0),shell=True,
                stdout=log_file,stderr=log_file))

    # Compile tex document with figures in it. 
    # Run twice to get references correct.
    if 'compile' in sys.argv[1:]:
        os.system("/usr/texbin/pdflatex gisr.tex")
        os.system("/usr/texbin/pdflatex gisr.tex")
