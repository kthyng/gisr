"""
Attempt at having a control-all file for oil modeling paper.
ADD EXAMPLE CALL
"""

import matplotlib
matplotlib.use("Agg") # set matplotlib to use the backend that does not require a windowing system
import numpy as np
import sys
import os
import op
import tracmass
import netCDF4 as netCDF
from mpl_toolkits.basemap import Basemap
import pdb
from matplotlib import delaunay
import matplotlib.pyplot as plt
import glob
from datetime import datetime, timedelta
from mpl_toolkits.basemap import Basemap
import time
import tracpy
import init
from scipy import ndimage
import projects
import matplotlib as mpl
import subprocess

# Units for time conversion with netCDF.num2date and .date2num
units = 'seconds since 1970-01-01'

if __name__ == "__main__":
    # Make sure necessary directories exist
    if not os.path.exists('logs'):
        os.makedirs('logs')
    process_queue = []
    poll_interval = 5.0 #600.0 # number of seconds to wait 
    arg_switches = {"sensitivity":False, "outer_f":False, "galv_b":False,
                    "bara_b":False, "dwh_f":False}#, "compile":False}
    # Check to see if any input arguments were on the command line
    if len(sys.argv) > 1:
        for arg in sys.argv[1:]:
            # Select out all of the flags sent in
            # Flags should be indicated using '--' in front
            key = arg[2:]
            # Set all input keys to true
            if key in arg_switches.keys():
                arg_switches[key] = True
            else:
                raise ValueError("Unknown test %s requested." % key)

    # Find number of cores to use
    # Use the environmental variable if set
    if os.environ.has_key('OMP_NUM_THREADS'):
        max_processes = int(os.environ['OMP_NUM_THREADS'])
    # Otherwise try to find out
    else:
        if 'Linux' in os.uname(): # Linux
            # Find number of cores on a Linux box
            temp = os.popen('less /proc/cpuinfo | grep processor').read()
            max_processes = temp.count('processor')

        elif 'Darwin' in os.uname(): # Mac
            # Find number of cores on a mac
            max_processes = os.popen('system_profiler | grep "Cores"').read()[-2]

    # Make list of process number
    proc_count = list(np.arange(max_processes))

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
    # The call assigns the process to a specific core by number
    # using `taskset [number]` and redirects the screen output
    # to a log file, and the process will run in the background
    for (key,value) in arg_switches.iteritems():
        if value:
            # Set up log file
            log_name.append('%s-%s' % (date,key))

            # Set up command
            cmd_list.append('python projects/%s.py' % key)
                            
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
            # pdb.set_trace()
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
