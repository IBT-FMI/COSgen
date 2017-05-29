import random

class MissingFunction(Exception):
	pass

def ga(population,functions,generations=100):
	
	if not hasattr(functions,'mutate'):
		raise MissingFunction("No 'mutate' function in 'functions'.")
	if not hasattr(functions,'cross_over'):
		raise MissingFunction("No 'cross_over' function in 'functions'.")
	if len(functions.fitness_measures) == 0:
		raise MissingFunction("No fitness measures in 'functions'.")
	population_size = len(population)
	n =  3 #int(population_size/2)	#number of sequences surviving each generation
	print('n = ',n)
	for seq in population:
		seq.fitness = functions.evaluate_fitness(seq)
	for i in range(generations):
		best_seqs = functions.find_best(population,n)
		population = best_seqs
		for i in range(population_size-n):
			population.append(functions.cross_over(best_seqs[random.randrange(n)],best_seqs[random.randrange(n)],mutation_probability=0.5))
			population[len(population)-1].fitness = functions.evaluate_fitness(population[len(population)-1])
	return population
