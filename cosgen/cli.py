#!/usr/bin/python3

try:
	from cosgen.function_crate import FunctionCrate
	from cosgen.algorithms import ga
	import cosgen.fitness_measures as fitness_measures
	from cosgen.sequence import Sequence
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

def cli_algorithm(population_size=20, library_size=20, storage_path='~/.cosgen', seqlength=100, nstimtypes=1, generations=10000, survivors=5, nimmigrants=4, hrflength=30, TR=1):

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
		f.write('TR = '+str(TR))

	fcts = FunctionCrate()
#	def design_mat(x):
#		return np.diag(x.l)
#	def cov_mat(x):
#		return np.identity(len(x[0]))
#	model = models.Model(design_mat, cov_mat)
	gamma_hrf = models.get_gamma_hrf(TR,hrflength,5,1,15,1,6)
	ar1_cov = models.get_ar1_cov(seqlength,0.5)
#	import math
#	gamma_hrf = np.zeros(40)
#	gamma_hrf[0:len(gamma_hrf):2]=1
#	ar1_cov = np.identity(139)
	extra_evs = np.empty((seqlength,2))
	extra_evs[:,0]=np.ones(seqlength)
	extra_evs[:,1]=np.linspace(-0.5,0.5,seqlength)
	model = models.DetectionModel(gamma_hrf, err_cov_mat=ar1_cov, filterfunc=partial(models.gaussian_highpass,sigma=225),extra_evs=extra_evs)
#	basis_set = models.get_FIR_basis_set(hrflength)
#	model = models.EstimationModel(basis_set,err_cov_mat=ar1_cov)
	fcts.add_fitness_measure('cov',partial(fitness_measures.estimator_variance,model=model,optimality='a'))
	fcts.set_mutate(mutate)
	fcts.set_cross_over(cross_over)
	fcts.set_generate_immigrants(partial(generate_immigrants, seqlen=seqlength, nstimtypes=nstimtypes, block_size=10))
	statistics = Statistics(storage_path)
	population = [Sequence(seqlength,nstimtypes) for i in range(population_size-1)]
	population.append(Sequence(seqlength,nstimtypes,'block',block_size=10))
	population = ga(population,fcts,generations,survivors,nimmigrants,statistics)
	for i,seq in enumerate(fcts.find_best(population,library_size)):
		seq.dump(storage_path,index=i,TR=TR)
		print(seq.l, seq.fitness)

def main():
	argh.dispatch_command(cli_algorithm)

if __name__ == '__main__':
	main()
