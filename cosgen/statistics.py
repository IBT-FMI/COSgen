from matplotlib import pyplot as plt

class Statistics:
	def __init__(self):
		self.max_fitness = []
		self.average_fitness = []
		self.population_correlation = []

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
		# TODO add population correlation

	def show(self):
		plt.plot(self.max_fitness)
		plt.title('Max fitness')
		plt.show()
		plt.plot(self.average_fitness)
		plt.title('Average fitness')
		plt.show()
