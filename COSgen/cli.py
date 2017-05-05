#!/usr/bin/python3

try:
	from COSgen.functions import functions
	from COSgen.algorithm import ga
	import COSgen.fitness_measures as fitness_measures
	from COSgen.sequence import sequence
	from COSgen.mutate_functions import mutate
	from COSgen.cross_over_functions import cross_over
except ImportError:
	from functions import functions
	from algorithm import ga
	import fitness_measures
	from sequence import sequence
	from mutate_functions import mutate
	from cross_over_functions import cross_over

import argh


def cli_algorithm(population_size=200, library_size=10, storage_path='~/.COSgen/sequences'):
	fcts = functions()
	fcts.add_fitness_measure('test',fitness_measures.test)
	fcts.set_mutate(mutate)
	fcts.set_cross_over(cross_over)
	population = [sequence() for i in range(population_size)]
	population = ga(population,fcts,generations=30)
	for seq in fcts.find_best(population,library_size):
		seq.dump(storage_path)

def main():
	argh.dispatch_command(cli_algorithm)

if __name__ == '__main__':
	main()
