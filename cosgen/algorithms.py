"""This file contains the genetic algorithm. Any extra algorithms should
be added in this file. """

import numpy as np
import math

class MissingFunction(Exception):
	"""
	Error raised if the function crate passed to 'ga' as 'functions'
	does not contain all necessary functions for the execution of the
	genetic algorithm.
	"""
	pass

def ga(population,functions,generations,nsurvive,nimmigrants,stat,print_freq=100):
	"""
	Run genetic algorithm.

	This function runs a genetic algorithm on `population` with the
	given arguments.

	Parameters
	---------
	population : list of cosgen.sequence.Sequence objects
	    Initial population.
	functions : cosgen.function_crate.FunctionCrate object
	    This object has to have at least a mutate and cross_over
	    function as well as one fitness measure.
	generations : int
	    Number of generations(iterations) of the genetic algorithm.
	nsurvive : int
	    Number of survivors after each generation.
	nimmigrants : int
	    Number of immigrants in each generation.
	stat : cosgen.statistics.Statistics object
	    Logs properties of population over generations.
	print_freq : int
	    Print generation number every print_freq generations.

	Returns
	-------
	population : list of cosgen.sequence.Sequence objects
	    Population after genetic algorithm.
	"""
	if not hasattr(functions,'mutate'):
		raise MissingFunction("No 'mutate' function in 'functions'.")
	if not hasattr(functions,'cross_over'):
		raise MissingFunction("No 'cross_over' function in 'functions'.")
	if len(functions.fitness_measures) == 0:
		raise MissingFunction("No fitness measures in 'functions'.")

	population_size = len(population)
	for seq in population:
		seq.fitness = functions.evaluate_fitness(seq)
	for gen in range(generations):
		if gen%print_freq == 0:
			print('Generation '+str(gen)+' of '+str(generations))
		stat.add(population)
		best_seqs = functions.find_best(population,nsurvive)
		population = best_seqs
		p = np.array([s.fitness for s in best_seqs])
		p = p/p.sum() #probability for choosing sequence as parent
		idxs = range(nsurvive) #indices for parent selection
		for i in range(population_size-nsurvive):
			parents = np.random.choice(idxs, 2, False, p)
			population.append(
				functions.mutate(
					functions.cross_over(best_seqs[parents[0]],best_seqs[parents[1]]),
					0.1)
			)
		population.extend(functions.generate_immigrants(nimmigrants))
		for seq in population:
			if math.isnan(seq.fitness):
				seq.fitness = functions.evaluate_fitness(seq)
	stat.gen_plot()
	return population
