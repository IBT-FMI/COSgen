import random
try:
	from cosgen.statistics import Statistics
except ImportError:
	from statistics import Statistics

class MissingFunction(Exception):
	pass

def ga(population,functions,generations,nsurvive):
	
	if not hasattr(functions,'mutate'):
		raise MissingFunction("No 'mutate' function in 'functions'.")
	if not hasattr(functions,'cross_over'):
		raise MissingFunction("No 'cross_over' function in 'functions'.")
	if len(functions.fitness_measures) == 0:
		raise MissingFunction("No fitness measures in 'functions'.")

	stat = Statistics()
	population_size = len(population)
	for seq in population:
		seq.fitness = functions.evaluate_fitness(seq)
	for i in range(generations):
		stat.add(population)
		best_seqs = functions.find_best(population,nsurvive)
		population = best_seqs
		for i in range(population_size-nsurvive):
			population.append(
				functions.mutate(
					functions.cross_over(best_seqs[random.randrange(nsurvive)],best_seqs[random.randrange(nsurvive)]),
					0.1)
			)
			population[len(population)-1].fitness = functions.evaluate_fitness(population[len(population)-1])
	stat.show()
	return population
