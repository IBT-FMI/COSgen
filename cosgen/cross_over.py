try:
	import cosgen.sequence as sequence
except ImportError:
	import sequence

import random

def cross_over(sequence1, sequence2, mutation_rate=0.3):
	length = min(len(sequence1.l),len(sequence2.l))
	r = random.randrange(length)
	l = sequence1.l[:r]
	l.extend(sequence2.l[r:])
	nstimtypes = max(sequence1.nstimtypes,sequence2.nstimtypes)
	return sequence.sequence(nstimtypes=nstimtypes,l=l)
