#!/usr/bin/python3

try:
	from cosgen.functions import functions
	from cosgen.algorithms import ga
	import cosgen.fitness_measures as fitness_measures
	from cosgen.sequence import sequence
	from cosgen.mutate import mutate
	from cosgen.cross_over import cross_over
	from cosgen.immigrants import generate_immigrants
	import cosgen.models
except ImportError:
	from functions import functions
	from algorithms import ga
	import fitness_measures
	from sequence import sequence
	from mutate import mutate
	from cross_over import cross_over
	from immigrants import generate_immigrants
	import cosgen.models

import argh
from os.path import expanduser
from functools import partial
import numpy as np

def cli_algorithm(population_size=200, library_size=10, storage_path='~/.cosgen/sequences'):
	storage_path = expanduser(storage_path)
	fcts = functions()
	def design_mat(x):
		return np.diag(x.l)
	def cov_mat(x):
		return np.identity(len(x[0]))
	model = cosgen.models.Model(design_mat, cov_mat)
	fcts.add_fitness_measure('test',partial(fitness_measures.estimator_variance,model=model,optimality='a'))
	fcts.set_mutate(mutate)
	fcts.set_cross_over(cross_over)
	fcts.set_generate_immigrants(generate_immigrants)
	population = [sequence(30,1) for i in range(population_size)]
	population = ga(population,fcts,generations=30)
	for seq in fcts.find_best(population,library_size):
		seq.dump(storage_path)

def main():
	argh.dispatch_command(cli_algorithm)

if __name__ == '__main__':
	main()
