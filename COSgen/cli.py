#!/usr/bin/python3

try:
	from COSgen.functions import functions
	from COSgen.algorithms import ga
	import COSgen.fitness_measures as fitness_measures
	from COSgen.sequence import sequence
	from COSgen.mutate import mutate
	from COSgen.cross_over import cross_over
	from COSgen.immigrants import generate_immigrants
except ImportError:
	from functions import functions
	from algorithms import ga
	import fitness_measures
	from sequence import sequence
	from mutate import mutate
	from cross_over import cross_over
	from immigrants import generate_immigrants

import argh
from os.path import expanduser

def cli_algorithm(population_size=200, library_size=10, storage_path='~/.COSgen/sequences'):
	storage_path = expanduser(storage_path)
	fcts = functions()
	fcts.add_fitness_measure('test',fitness_measures.test)
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
