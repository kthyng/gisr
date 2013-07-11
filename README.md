#GISR Project: TX-LA shelf drifter tracking

Run the GISR project and save the standard output to a log file at the command prompt with `python2.7 run.py > log.txt`

Use flags in run.py to turn on/off various pieces of the work.

Even better, multiple instances of the simulation can be run on different cores if you have a multi-core Linux machine. To do this, use `taskset`:
`taskset 3 python2.7 run.py > temp_log.txt &`
This command runs the run.py script for a specific project in the background (due to the &), redirects the output from the screen to the text file 'temp_log.txt', and runs the process on core 3. This command can then be used for other instances by choosing a different core after the command `taskset`.