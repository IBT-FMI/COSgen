from functools import partial

class WrongOrderError(Exception):
	pass

class RmAttrError(Exception):
	pass

class OverwriteAttrError(Exception):
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
		if not name in self.fitness_measure:
			print('Cannot remove fitness measure "{0}" because it does not exist.')
			return
		del self.fitness_measures[name]

	def set_mutate(self,function):
		if hasattr('mutate'):
			raise OverwriteAttrError('Use del_mutate before setting a new mutate')
		setattr(self,'mutate',function)

	def del_mutate(self):
		if not hasattr('mutate'):
			raise RmAttrError('Cannot remove attribute that was not set.') 
		delattr(self,'mutate')

	def set_cross_over(self,function):
		if hasattr('cross_over'):
			raise OverwriteAttrError('Use del_cross_over before setting a new cross over')
		setattr(self,'cross_over',function)

	def del_cross_over(self):
		if not hasattr('cross_over'):
			raise RmAttrError('Cannot remove attribute that was not set.') 
		delattr(self,'cross_over')

	def set_generate_immigrants(self,function):
		if hasattr('generate_immigrants'):
			raise OverwriteAttrError('Use del_generate_immigrants before setting a new generate immigrants')
		if not hasattr(self,'cross_over'):
			raise WrongOrderError('cross_over function has to be set before set_generate_imigrants can be used')
		setattr(self,'generate_immigrants',partial(function,cross_over_fct=self.cross_over))

	def del_generate_immigrants(self):
		if not hasattr('generate_immigrants'):
			raise RmAttrError('Cannot remove attribute that was not set.') 
		delattr(self,'generate_immigrants')

	def set_model(self,model):
		if hasattr('model'):
			raise OverwriteAttrError('Use del_model before setting a new model')
		setattr(self,'model',model)

	def del_model(self)
		if not hasattr('model'):
			raise RmAttrError('Cannot remove attribute that was not set.') 
		delattr(self,'model')
