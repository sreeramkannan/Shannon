import multiprocessing
import os
import run_MB_SF_fn

def run_process(cmd):
    os.system(cmd)

def run_cmds(cmds,noJobs):
	#cmds is a tuple of strings
	#noJobs is a an integer with the no of Jobs to run in parallel
	p = multiprocessing.Pool(noJobs)
	p.map(run_process,cmds)

def run_MB_SF(cmds):
	p = multiprocessing.Pool(noJobs)
	p.map(run_MB_SF_fn.run_MB_SF_fn,cmds)