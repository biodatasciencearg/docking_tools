#Based on the original script from https://gist.github.com/joezuntz/7c590a8652a1da6dc4c9 and modified by Elias Lopez 
# for the Autodock implementation   
# Contact: msn.eliaslopez@gmail.com
#####################################################################################################################
#Suppose you have a collection of tasks, which in this example I'll assume is just running a function f.
#If these tasks are completely separate and independent the most then you can parallelize them easily.
#In this gist I'll show the simplest possible way to do this using mpi4py.
#There are better ways to do this, in particular if the tasks vary significantly in time taken to run.
####################################################################################################################
#!/usr/bin/env python


import os, sys
from optparse import OptionParser
from optparse import Option, OptionValueError
import subprocess
import mpi4py.MPI
import glob
from os.path import basename, splitext
#Funcion trivial para paralelizar
def f(i):
  #write some info to it
  #print("Esto andara?? infile %s \n"%(i))
  var=i
  filename =  glob.glob(var + "/*.dpf")[0]
  filename2 = os.path.splitext(basename(filename))[0]
  return_code = subprocess.call("cd " + var + "&& autodock4 -p "+filename2+".dpf -l " + filename2 + ".dlg.pdb", shell=True) 



VERSION = '1.0'
# Simplemente defino clases para parsear la entrada.
class MultipleOption(Option):
    ACTIONS = Option.ACTIONS + ("extend",)
    STORE_ACTIONS = Option.STORE_ACTIONS + ("extend",)
    TYPED_ACTIONS = Option.TYPED_ACTIONS + ("extend",)
    ALWAYS_TYPED_ACTIONS = Option.ALWAYS_TYPED_ACTIONS + ("extend",)

    def take_action(self, action, dest, opt, value, values, parser):
        if action == "extend":
            values.ensure_value(dest, []).append(value)
        else:
            Option.take_action(self, action, dest, opt, value, values, parser)


def main():
    PROG = os.path.basename(os.path.splitext(__file__)[0])
    long_commands = ('categories')
    short_commands = {'cat':'categories'}
    description = 'Description of command...\n -f     folders\n							Developed by Elias Lopez (2016) contact:msn.eliaslopez@gmail.com'
    parser = OptionParser(option_class=MultipleOption,
                          usage='usage: %prog [OPTIONS] COMMAND [folders]',
                          version='%s %s' % (PROG, VERSION),
                          description=description)
    parser.add_option('-f', '--folders', 
                      action="extend", type="string",
                      dest='folders', 
                      metavar='FOLDERS', 
                      help='')

    if len(sys.argv) == 1:
        parser.parse_args(['--help'])

    OPTIONS, args = parser.parse_args()
    #A list of all the tasks to do.  In your case you will probably build this task list in a more complex way.
    #You don't even need to build it in advance for this approach to work
    d = str(OPTIONS).split(":")[1].strip("}, ], [, ', '")
    task_list = args   
    task_list.append(str(d))
    #And now moving on the parallel version
 
    #mpi4py has the notion of a "communicator" - a collection of processors
    #all operating together, usually on the same program.  Each processor 
    #in the communicator is identified by a number, its rank,  We'll use that
    #number to split the tasks
    #find out which number processor this particular instance is,
    #and how many there are in total
    rank = mpi4py.MPI.COMM_WORLD.Get_rank()
    size = mpi4py.MPI.COMM_WORLD.Get_size()
    name = mpi4py.MPI.Get_processor_name()
    #parallelized version
    #the enumerate function gives us a number i in addition
    #to the task.  (In this specific case i is the same as task!  But that's
    #not true usually)
    for i,task in enumerate(task_list):
    	#This is how we split up the jobs.
	#The % sign is a modulus, and the "continue" means
	#"skip the rest of this bit and go to the next time
	#through the loop"
	# If we had e.g. 4 processors, this would mean
	# that proc zero did tasks 0, 4, 8, 12, 16, ...
	# and proc one did tasks 1, 5, 9, 13, 17, ...
	# and do on.
	if i%size!=rank: continue
	print "Task number %d (%s) is running by processor %d of %d on %s" % (i, task, rank, size, name)
    	f(task)
        print "Task number %d (%s) being done by processor %d of %d on %s" % (i, task, rank, size, name)


if __name__ == '__main__':
    main()
