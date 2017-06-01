from functools import partial
import random

class WrongOrderError(Exception):
	pass

class RmAttrError(Exception):
	pass

class OverwriteAttrError(Exception):
	pass

def partition(population, left, right, pivotIndex):
	"""code from https://rosettacode.org/wiki/Quickselect_algorithm#Python """
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
	"""code from https://rosettacode.org/wiki/Quickselect_algorithm#Python """
	""""Returns the k-th smallest, (k >= 0), element of population within population[left:right+1] inclusive."""
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
		self.fitness_measures = {}

	def evaluate_fitness(self,sequence):
		return sum([self.fitness_measures[f](sequence) for f in self.fitness_measures])

	@staticmethod
	def find_best(population,n):
		nidx = quickselect(population,0,len(population)-1,len(population)-n)
		return population[len(population)-n:]

	def add_fitness_measure(self,name,function):
		if name in self.fitness_measures:
			print("A fitness measure called {0} already "
			      "exists. It will be overwritten.".format(name))
		self.fitness_measures[name] = function

	def remove_fitness_measure(self,name):
		if not name in self.fitness_measure:
			print('Cannot remove fitness measure "{0}" because it does not exist.'.format(name))
			return
		del self.fitness_measures[name]

	def set_mutate(self,function):
		if hasattr(self,'mutate'):
			raise OverwriteAttrError('Use del_mutate before setting a new mutate')
		setattr(self,'mutate',function)

	def del_mutate(self):
		if not hasattr(self,'mutate'):
			raise RmAttrError('Cannot remove attribute that was not set.') 
		delattr(self,'mutate')

	def set_cross_over(self,function):
		if hasattr(self,'cross_over'):
			raise OverwriteAttrError('Use del_cross_over before setting a new cross over')
		setattr(self,'cross_over',function)

	def del_cross_over(self):
		if not hasattr(self,'cross_over'):
			raise RmAttrError('Cannot remove attribute that was not set.') 
		delattr(self,'cross_over')

	def set_generate_immigrants(self,function):
		if hasattr(self,'generate_immigrants'):
			raise OverwriteAttrError('Use del_generate_immigrants before setting a new generate immigrants')
		if not hasattr(self,'cross_over'):
			raise WrongOrderError('cross_over function has to be set before set_generate_imigrants can be used')
		setattr(self,'generate_immigrants',partial(function,cross_over_fct=self.cross_over))

	def del_generate_immigrants(self):
		if not hasattr(self,'generate_immigrants'):
			raise RmAttrError('Cannot remove attribute that was not set.') 
		delattr(self,'generate_immigrants')
