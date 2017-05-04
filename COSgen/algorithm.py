def find_best(population):
	pass

class MissingFunction(Exception):
	pass

def ga(population,functions,generations=100):
	
	if not hasattr(functions,'mutate'):
		raise MissingFunction("No 'mutate' function in 'functions'.")
	if not hasattr(functions,'cross_over'):
		raise MissingFunction("No 'cross_over' function in 'functions'.")
	if len(functions.fitness_measures) == 0:
		raise MissingFunction("No fitness measures in 'functions'.")

	for seq in population:
		seq.fitness = functions.evaluate_fitness(seq)
	for i in range(generations):
		best_seqs = find_best(population)
		#cross over, mutate and update the fitnesses
