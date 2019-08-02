#!/usr/bin/env python

"""
Parallel implementation of gadrad using mpi4py. This delegates tasks to
a number of cores thus allowing the genetic algorithm to find solutions much faster
"""

from mpi4py import MPI
import numpy as np
import radmc3dpy
import GA
from datetime import datetime
import time
import os
import shutil
import params
import math
import radmc3dpy
import traceback

pop_size = params.ga_params()['pop_size']
max_gen = params.ga_params()['max_gen']
ga_num_max = params.ga_params()['ga_num']
local_output = params.files()['outputs']
global_input = params.files()['inputs']
fit_lim = params.ga_params()['fit_lim']

def enum(*sequential, **named):
    """Handy way to fake an enumerated type in Python
    http://stackoverflow.com/questions/36932/how-can-i-represent-an-enum-in-python
    """
    enums = dict(zip(sequential, range(len(sequential))), **named)
    return type('Enum', (), enums)

# Define MPI message tags
tags = enum('READY', 'DONE', 'EXIT', 'START', 'GEN')

# Initializations and preliminaries
comm = MPI.COMM_WORLD   # get MPI communicator object
size = comm.size        # total number of processes
rank = comm.rank        # rank of this process
status = MPI.Status()   # get MPI status object


gen = 0
numga = 0
#bestfit=float('inf')

child_chromes = [None]*(pop_size+1)


if rank == 0:
    t0 = time.time()
    import params

while numga <= ga_num_max:
    bestfit=float('inf')
    numga += 1
    gen = 0


    if rank == 0:
        numga_str = str(datetime.utcnow().strftime('%y%m%d_%H%M%S%f'))
        bestfit=float('inf')
    else:
        numga_str = None
        bestfit=None
    numga_str = comm.bcast(numga_str,root=0)

    if rank != 0:
        dir_rank = str(rank).zfill(6)
	shutil.copytree(global_input,local_output+'/'+numga_str+'/'+dir_rank,symlinks=False)
	os.chdir(local_output+'/'+numga_str+'/'+dir_rank)


    while (gen <= max_gen):
        bestfit = comm.bcast(bestfit,root=0)
        if bestfit < fit_lim:
            break

	gen += 1
	if rank == 0:
            # Master process executes code below
	    #import params
	    t1 = time.time()
            tasks = range(pop_size)
	    task_index = 0
	    num_workers = size - 1
	    closed_workers = 0
            fit_chromes = [None]*(len(params.var_params())*2+2)

            while closed_workers < num_workers:
		data = comm.recv(source=MPI.ANY_SOURCE, tag=MPI.ANY_TAG, status=status)
		source = status.Get_source()
		tag = status.Get_tag()
		if tag == tags.READY:
			# Worker is ready, so send it a task
	            if task_index < len(tasks):
	                fdata = [child_chromes[task_index]]
			#outfile = open('fdata','w')
                        #print >> outfile, fdata
                        comm.send(fdata,dest=source,tag=tags.START)
			task_index += 1
		    else:
			comm.send(None, dest=source, tag=tags.EXIT)
		elif tag == tags.DONE:
		    #print data
                    fit_chromes = np.vstack((fit_chromes,data))
		elif tag == tags.EXIT:
		    closed_workers += 1

	    fit_chromes = np.delete(fit_chromes,0,0)
	    bestfit = min(fit_chromes[:,0])
	    child_chromes = GA.children(fit_chromes)
            #outfile = open('child_chromes','w')
            #print >> outfile, child_chromes

	else:
	    while True:
		comm.send(None, dest=0, tag=tags.READY)
		task = comm.recv(source=0, tag=MPI.ANY_TAG, status=status)
		tag = status.Get_tag()
		chromo_keys = GA.chromo_keys()
                if tag == tags.START:
			# Do the work here
		    if gen == 1:
			chromo_init   = GA.init()
                        chromo_values = chromo_init[0]

			chromo = dict(zip(chromo_keys,chromo_values))

                        try:
                            fitness = [radmc3dpy.runtask(gen,rank,chromo)]

                        except Exception:
			    fitness = [float("inf")]
			    print "consider changing parameter range"
                            s = traceback.format_exc()
                            print s

                        fitness.extend([gen])
                        fitness.extend(chromo_values)
                        fitness.extend(chromo_init[1])

			comm.send(fitness, dest=0, tag=tags.DONE)

		    else:

                        chromo_values = task[0]
                        #test = task[1]
			chromo = dict(zip(chromo_keys,chromo_values))

                        try:
	                    fitness = [radmc3dpy.runtask(gen,rank,chromo)]
			except:
			    fitness = [float("inf")]
			    print "consider changing parameter range"

                        fitness.extend([gen])
                        fitness.extend(chromo_values)

			comm.send(fitness, dest=0, tag=tags.DONE)

		elif tag == tags.EXIT:
	            break

	    comm.send(None, dest=0, tag=tags.EXIT)

	comm.Barrier()
