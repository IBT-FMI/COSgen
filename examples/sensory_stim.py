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
import pandas as pd
from functools import partial

storage_path = '~/.cosgen'
seqlength = 1490
TR = 1
nstimtypes = 1
population_size = 20

ecm = models.get_ar1_cov(seqlength,0.78)	#ecm = estimator covariance matrix
whitening_mat = np.linalg.inv(np.linalg.cholesky(ecm))

#def design_matrix(sequence):
#	hrf = np.load('/home/wguest/hrf.npy')
#	hrflength = len(hrf)
#	ls = len(sequence.l)
#	DM = np.empty((ls,4))
#	DM[:,0]=np.ones(seqlength)
#	DM[:,1]=np.linspace(-0.5,0.5,seqlength)
#	X = np.array([sequence.l == i for i in range(1,sequence.nstimtypes+1)],dtype=int)
#	blocks = sequence.get_block_representation()
#	Xconfoundregressor = np.zeros((1,sequence.seqlen))
#	for idx, val in enumerate(blocks):
#		Xconfoundregressor[0,val[0]:val[1]] = len(blocks)-idx
#	DM[:,2:3] = np.transpose(np.apply_along_axis(lambda m: models.orthogonalize(DM[:,0:2],models.gaussian_highpass(np.convolve(m,hrf)[0:ls])), axis=1, arr=X))
#	DM[:,3:] = np.transpose(np.apply_along_axis(lambda m: models.orthogonalize(DM[:,0:3],models.gaussian_highpass(np.convolve(m,hrf)[0:ls])), axis=1, arr=Xconfoundregressor))
#	return DM
#
#def cov_beta(X):
#	Z = whitening_mat*X
#	try:
#		Zpinv = np.linalg.pinv(Z)
#	except np.linalg.linalg.LinAlgError:
#		print('X:')
#		print(X)
#		print('K:')
#		print(self.whitening_mat)
#		print('Z:')
#		print(Z)
#	return Zpinv * np.transpose(Zpinv)
def main():
	fsl_basis_set_file = '/usr/share/fsl/etc/default_flobs.flobs/hrfbasisfns.txt'
	hrf = np.loadtxt(fsl_basis_set_file, ndmin=2)
	#hrf = models.get_FIR_basis_set(3)
	model = models.EstimationModel(hrf, whitening_mat=whitening_mat)

	fc = FunctionCrate()
	c = np.zeros(6)
	c[0]=1
	c[1]=1
	c[2]=1
	c[3]=-1
	c[4]=-1
	c[5]=-1
	contrast = np.matrix(c)
	fc.add_fitness_measure('est_var', partial(fitness_measures.estimator_variance, model=model, optimality='d', contrast=contrast))
	blocksize = estimate_optimal_block_size(seqlength, fc)
	print('blocksize', blocksize)
	return
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
