#!/usr/bin/python3

try:
	from cosgen.function_crate import FunctionCrate
	from cosgen.algorithms import ga
	import cosgen.fitness_measures as fitness_measures
	from cosgen.sequence import Sequence, estimate_optimal_block_size
	from cosgen.mutate import mutate
	from cosgen.cross_over import cross_over
	from cosgen.immigrants import generate_immigrants
	import cosgen.models as models
	from cosgen.statistics import Statistics
except ImportError:
	from function_crate import FunctionCrate
	from algorithms import ga
	import fitness_measures
	from sequence import Sequence
	from mutate import mutate
	from cross_over import cross_over
	from immigrants import generate_immigrants
	import models
	from statistics import Statistics

import argh
import os.path
import os
from functools import partial
import datetime
import numpy as np

def cli_algorithm(population_size=20, library_size=20, storage_path='~/.cosgen', seqlength=100, nstimtypes=1, generations=10000, survivors=5, nimmigrants=4, hrflength=30, TR=1, model_type='detection', autoregression=0.5, baseline='auto'):
	"""
	Run default optimization.

	This function attempts to optimize find a good sequence with only few input parameters.
	It uses a generic design matrix construction and a first order autoregressive model.

	Parameters
	----------
	population_size : int
	    Size of the population used for the genetic algorithm.
	library_size : int
	    Number of sequences included in the output. The library_size best
	    sequences are saved.
	storage_path : string
	    Path to folder where output is saved.
	seqlength : int
	    Length of the sequence that is optimized.
	nstimtypes : int
	    Number of possible stimulus types in the sequence.
	generations : int
	    Number of iterations the genetic algorithm performs.
	survivors : int
	    Number of sequences that are carried over from the previous iteration.
	nimmigrants : int
	    Number of new (randomly generated) sequences add in each iteration.
	hrflength : int
	    Length of the HRF.
	TR : float
	    Repetition time of scans in seconds.
	model_type : string
	    Can be 'detection' in order to optimize for contrast detection
	    or 'estimation' for HRF estimation.
	autoregression : float
	    Coefficient of the fist order autoregressive noise model.
	baseline : string or int
	    Number of baseline TR before the stimulation starts. Can be 'auto'.
	"""
	storage_path = os.path.expanduser(storage_path)
	storage_path = os.path.join(storage_path,'{:%Y%m%d%H%M%S}'.format(datetime.datetime.now()))
	try:
		os.makedirs(storage_path)
	except OSError as exc:
		if exc.errno == errno.EEXIST and os.path.isdir(path):
			pass
		else:
			raise
	with open(os.path.join(storage_path,'parameters.txt'),'w+') as f:
		f.write('Population size = '+str(population_size)+'\n')
		f.write('Library size = '+str(library_size)+'\n')
		f.write('Sequence length = '+str(seqlength)+'\n')
		f.write('Number of stimulus types = '+str(nstimtypes)+'\n')
		f.write('Number of generations = '+str(generations)+'\n')
		f.write('Number of survivors = '+str(survivors)+'\n')
		f.write('Number of immigrants = '+str(nimmigrants)+'\n')
		f.write('HRF length = '+str(hrflength)+'\n')
		f.write('TR = '+str(TR)+'\n')
		f.write('Autoregression = '+str(autoregression)+'\n')
		f.write('Baseline = '+str(baseline))

	fcts = FunctionCrate()
	gamma_hrf = models.get_gamma_hrf(TR,hrflength)
	ar1_cov = models.get_ar1_cov(seqlength,autoregression)
	extra_evs = np.empty((seqlength,2))
	extra_evs[:,0]=np.ones(seqlength)
	extra_evs[:,1]=np.linspace(-0.5,0.5,seqlength)
	if model_type == 'detection':
		model = models.DetectionModel(gamma_hrf, err_cov_mat=ar1_cov, filterfunc=partial(models.gaussian_highpass,sigma=225),extra_evs=extra_evs)
	elif model_type == 'estimation':
		basis_set = models.get_FIR_basis_set(hrflength)
		model = models.EstimationModel(basis_set,err_cov_mat=ar1_cov)
	fcts.add_fitness_measure('cov',partial(fitness_measures.estimator_variance,model=model,optimality='a'))
	fcts.set_mutate(mutate)
	fcts.set_cross_over(cross_over)
	optimal_block_size = estimate_optimal_block_size(seqlength, fcts)
	if baseline=='auto':
		baseline = 2*optimal_block_size
	fcts.set_generate_immigrants(partial(generate_immigrants, seqlen=seqlength, nstimtypes=nstimtypes, block_size=optimal_block_size))
	statistics = Statistics(storage_path)
	population = [Sequence(seqlength,nstimtypes) for i in range(population_size-1)]
	population.append(Sequence(seqlength,nstimtypes,'block',block_size=10))
	population = ga(population,fcts,generations,survivors,nimmigrants,statistics)
	for i,seq in enumerate(fcts.find_best(population,library_size)):
		seq.add_baseline(baseline)
		seq.dump(storage_path,index=i,TR=TR)
		print(seq.l, seq.fitness)

def main():
	argh.dispatch_command(cli_algorithm)

if __name__ == '__main__':
	main()
