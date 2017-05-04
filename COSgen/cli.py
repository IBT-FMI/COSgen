#!/usr/bin/python3

from functions import functions
from algorithm import ga
import fitness_measures
from sequence import sequence
from mutate_functions import mutate
from cross_over_functions import cross_over

if __name__ == "__main__":

	fcts = functions()
	fcts.add_fitness_measure('test',fitness_measures.test)
	fcts.set_mutate(mutate)
	fcts.set_cross_over(cross_over)
	population = [sequence() for i in range(100)]

	ga(population,fcts,generations=30)
