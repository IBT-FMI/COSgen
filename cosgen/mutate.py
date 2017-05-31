import random

def mutate(sequence,mutation_prob):
	length = len(sequence.l)
	idxs = random.sample(range(length),int(length*mutation_prob))
	for i in idxs:
		sequence.l[i] = random.randrange(sequence.nstimtypes+1)
	return sequence
