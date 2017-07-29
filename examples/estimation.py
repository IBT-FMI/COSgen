#!/usr/bin/python3

from cosgen.function_crate import FunctionCrate
from cosgen.algorithms import ga
import cosgen.fitness_measures as fitness_measures
from cosgen.sequence import Sequence, estimate_optimal_block_size
from cosgen.mutate import mutate
from cosgen.cross_over import cross_over
from cosgen.immigrants import generate_immigrants
import cosgen.models as models
from cosgen.statistics import Statistics

import numpy as np
from functools import partial
from scipy.ndimage.filters import gaussian_filter1d

storage_path = '~/.cosgen/estimation_FIR10'
seqlength = 1490
TR = 1
nstimtypes = 1
population_size = 20

ecm = models.get_ar1_cov(seqlength,0.9398)	#ecm = estimator covariance matrix
whitening_mat = np.linalg.inv(np.linalg.cholesky(ecm))

filter_func = lambda x: gaussian_filter1d(models.gaussian_highpass(x),2)
extra_evs = np.empty((seqlength,2))
extra_evs[:,0]=np.ones(seqlength)
extra_evs[:,1]=np.linspace(-0.5,0.5,seqlength)

def main():
	n_basis_vectors= 10
	basis_set = models.get_FIR_basis_set(n_basis_vectors)
	model = models.EstimationModel(basis_set, whitening_mat=whitening_mat, filterfunc=filter_func, nonlincorrection=models.squashing_function, extra_evs=extra_evs)

	fc = FunctionCrate()
	c = np.zeros((n_basis_vectors,n_basis_vectors+2)) #+2 for extra explanatory variables
	for i in range(n_basis_vectors):
		c[i,i+2]=1
	c=np.matrix(c)
	print(c)
	fc.add_fitness_measure('est_var', partial(fitness_measures.estimator_variance, model=model, optimality='a', contrast=c))
	blocksize = estimate_optimal_block_size(seqlength, fc)
	print('blocksize', blocksize)
	fc.set_mutate(mutate)
	fc.set_cross_over(cross_over)
	fc.set_generate_immigrants(partial(generate_immigrants, seqlen=seqlength, nstimtypes=nstimtypes, block_size=blocksize))
	statistics = Statistics(storage_path)
	population = [Sequence(seqlength,nstimtypes) for i in range(population_size-1)]
	population.append(Sequence(seqlength,nstimtypes,'block',block_size=blocksize))
	population = ga(population,fc,generations=10000,nsurvive=5,nimmigrants=4,stat=statistics)
	fc.find_best(population,1)[0].dump(storage_path,TR=TR)

if __name__ == "__main__":
	main()
