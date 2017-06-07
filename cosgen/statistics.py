from matplotlib import pyplot as plt
import numpy as np

class Statistics:
	def __init__(self):
		self.max_fitness = []
		self.average_fitness = []
		self.population_variance = []

	def add(self,population):
		maxfitness = 0
		avefitness = 0
		for s in population:
			if s.fitness > maxfitness:
				maxfitness = s.fitness
			avefitness += s.fitness
		avefitness = avefitness/len(population)
		self.max_fitness.append(maxfitness)
		self.average_fitness.append(avefitness)
		self.population_variance.append(1 - np.mean(np.absolute(np.corrcoef([s.l for s in population],rowvar=False))))

	def show(self):
		f, axarr = plt.subplots(3,sharex=True)
		axarr[0].plot(self.max_fitness)
		axarr[0].set_title('Max fitness')
		axarr[1].plot(self.average_fitness)
		axarr[1].set_title('Average fitness')
		axarr[2].plot(self.population_variance)
		axarr[2].set_title('Population variance')
		plt.show()
