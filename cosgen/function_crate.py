from functools import partial
import random
import inspect

class MissingAttrError(Exception):
	"""This error is raised if a FunctionCrates object misses an
	attribute to complete the requested operation."""
	pass

class WrongOrderError(Exception):
	"""This error is raised if functions are added to a
	FunctionCrate object in the wrong order."""
	pass

class RmAttrError(Exception):
	"""This error is raised if an attribute of a FunctionCrate object
	can not be removed because it does not exist."""
	pass

class OverwriteAttrError(Exception):
	"""This error is raised if an attribute of a FunctionCrate object
	already exist."""
	pass

def partition(population, left, right, pivotIndex):
	"""Helper function for quickselect.
	(Code from https://rosettacode.org/wiki/Quickselect_algorithm#Python)"""
	pivotValue = population[pivotIndex].fitness
	population[pivotIndex], population[right] = population[right], population[pivotIndex]  # Move pivot to end
	storeIndex = left
	for i in range(left, right):
	    if population[i].fitness < pivotValue:
	        population[storeIndex], population[i] = population[i], population[storeIndex]
	        storeIndex += 1
	population[right], population[storeIndex] = population[storeIndex], population[right]  # Move pivot to its final place
	return storeIndex

def quickselect(population, left, right, k):
	"""Returns the k-th smallest, (k >= 0), element of `population` within `population[left:right+1]` inclusive.
	Implementation of the quickselect algorithm from https://rosettacode.org/wiki/Quickselect_algorithm#Python."""
	while True:
		pivotIndex = random.randint(left, right)     # select pivotIndex between left and right
		pivotNewIndex = partition(population, left, right, pivotIndex)
		pivotDist = pivotNewIndex - left
		if pivotDist == k:
			return population[pivotNewIndex]
		elif k < pivotDist:
			right = pivotNewIndex - 1
		else:
			k -= pivotDist + 1
			left = pivotNewIndex + 1

class FunctionCrate:

	def __init__(self):
		"""This class holds all necessary functions for the genetic algorithm."""
		self.fitness_measures = {}

	def evaluate_fitness(self, sequence):
		"""
		Calculate overall fitness measure.

		This method calculates the sum of the return values of
		all fitness measure functions added to the instance of the class.

		Parameters
		----------
		sequence : cosgen.sequence.Sequence
		    Sequence for which the fitness is calculated.

		Returns
		-------
		float
		    Overall fitness.
		"""

		if len(self.fitness_measures) == 0:
			raise MissingAttrError("""No fitness measure in this FunctionCrate object {0}.
					Add fitness measures using the add_fitness_measure method.""".format(self))

		return sum([self.fitness_measures[f](sequence) for f in self.fitness_measures])

	@staticmethod
	def find_best(population, n):
		"""
		Find `n` best sequences in `population`.

		This method finds the `n` sequences with the highest fitness in the `population`.

		Parameters
		----------
		population : list of cosgen.sequence.Sequence
		    Population of sequences
		n : int
		    Number of sequences returned.

		Returns
		-------
		list of cosgen.sequence.Sequence
		    List of best sequences.
		"""
		nidx = quickselect(population, 0, len(population)-1, len(population)-n)
		return population[len(population)-n:]

	def add_fitness_measure(self, name, function):
		"""
		Add fitness measure function to object.

		This method adds a fitness measure function to the
		object, that is used in the `evaluate_fitness` method. The
		function must take a sequences as parameter and retrun a
		float.

		Parameters
		----------
		name : string
		    Name of the fitness measure.
		function : function
		    Fitness measure function.
		"""
		if name in self.fitness_measures:
			print("A fitness measure called {0} already "
			      "exists. It will be overwritten.".format(name))
		self.fitness_measures[name] = function

	def remove_fitness_measure(self, name):
		"""
		Remove fitness measure.

		This method removes a fitness measure previously added
		with the `add_fitness_measure` method.

		Parameters
		----------
		name : string
		    Name of the fitness measure to be removed.
		"""
		if not name in self.fitness_measure:
			print('Cannot remove fitness measure "{0}" because it does not exist.'.format(name))
			return
		del self.fitness_measures[name]

	def set_mutate(self, function):
		"""
		Set a mutate function.

		This method sets a mutate function used by the genetic
		algorithm. The function should normally take a
		`cosgen.sequence.Sequence` object as input and return a
		`cosgen.sequence.Sequence` object.

		Parameters
		----------
		function : function
		    Mutate function.
		"""
		if hasattr(self, 'mutate'):
			raise OverwriteAttrError('Use del_mutate before setting a new mutate')
		setattr(self, 'mutate', function)

	def del_mutate(self):
		"""
		Delete mutate function.

		This methods deletes the function added using the
		`set_mutate` method.
		"""
		if not hasattr(self, 'mutate'):
			raise RmAttrError('Cannot remove attribute that was not set.')
		delattr(self, 'mutate')

	def set_cross_over(self, function):
		"""
		Set a cross over function.

		This method sets a cross over function used by the
		genetic algorithm. The function should normally take two
		`cosgen.sequence.Sequence` objects as input and return a
		`cosgen.sequence.Sequence` object.

		Parameters
		----------
		function : function
		    Cross over function.
		"""
		if hasattr(self, 'cross_over'):
			raise OverwriteAttrError('Use del_cross_over before setting a new cross over')
		setattr(self,'cross_over', function)

	def del_cross_over(self):
		"""
		Delete cross over function.

		This methods deletes the function added using the
		`set_cross_over` method.
		"""
		if not hasattr(self, 'cross_over'):
			raise RmAttrError('Cannot remove attribute that was not set.')
		delattr(self, 'cross_over')

	def set_generate_immigrants(self,function):
		"""
		Set a generate immigrants function.

		This method sets a generate immigrants function used by
		the genetic algorithm. The function should normally return
		a list of `cosgen.sequence.Sequence` objects. `function`
		may uses a cross over function for the construction
		of immigrants. If `function` has a parameter
		'cross_over_fct' it is automatically set to the function
		previously set with `set_cross_over`.

		Parameters
		----------
		function : function
		    Generate immigrants function.
		"""
		if hasattr(self,'generate_immigrants'):
			raise OverwriteAttrError('Use del_generate_immigrants before setting a new generate immigrants')
		try:
			if 'cross_over_fct' in (inspect.getfullargspec(function).args + inspect.getfullargspec(function).kwonlyargs):
				if not hasattr(self, 'cross_over'):
					raise WrongOrderError('cross_over function has to be set before set_generate_imigrants can be used')
				setattr(self,'generate_immigrants', partial(function,cross_over_fct=self.cross_over))
			else:
				setattr(self,'generate_immigrants', function)
		except AttributeError:
			if 'cross_over_fct' in inspect.getargspec(function.func).args:
				if not hasattr(self, 'cross_over'):
					raise WrongOrderError('cross_over function has to be set before set_generate_imigrants can be used')
				setattr(self,'generate_immigrants', partial(function,cross_over_fct=self.cross_over))
			else:
				setattr(self,'generate_immigrants', function)

	def del_generate_immigrants(self):
		"""
		Delete generate immigrants function.

		This methods deletes the function added using the
		`set_generate_immigrants` method.
		"""
		if not hasattr(self,'generate_immigrants'):
			raise RmAttrError('Cannot remove attribute that was not set.')
		delattr(self,'generate_immigrants')
