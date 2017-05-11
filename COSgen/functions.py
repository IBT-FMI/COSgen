from functools import partial

class WrongOrderError(Exception):
	pass

class functions:

	def __init__(self):
		self.fitness_measures = {}

	def evaluate_fitness(self,sequence):
		return sum([self.fitness_measures[f](sequence) for f in self.fitness_measures])

	@staticmethod
	def find_best(population,n):
		#implement this properly!! sequence is not sorted!!
		return population[:n]

	def add_fitness_measure(self,name,function):
		if name in self.fitness_measures:
			print("A fitness measure called {0} already "
			      "exists. It will be overwritten.".format(name))
		self.fitness_measures[name] = function

	def remove_fitness_measure(self,name):
		del self.fitness_measures[name]

	def set_mutate(self,function):
		setattr(self,'mutate',function)

	def del_mutate(self,function):
		delattr(self,'mutate')

	def set_cross_over(self,function):
		setattr(self,'cross_over',function)

	def del_cross_over(self,function):
		delattr(self,'cross_over')

	def set_generate_immigrants(self,function):
		if not hasattr(self,'cross_over'):
			raise WrongOrderError('cross_over function has to be set before set_generate_imigrants can be used')
		setattr(self,'generate_immigrants',partial(function,cross_over_fct=self.cross_over))

	def del_generate_immigrants(self,function):
		delattr(self,'generate_immigrants')
