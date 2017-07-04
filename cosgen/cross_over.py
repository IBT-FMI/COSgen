"""This file contains the cross_over function. Any extra cross over
functions should be added in this file."""
try:
	import cosgen.sequence as sequence
except ImportError:
	import sequence

import random
import numpy as np

def cross_over(sequence1, sequence2):
	"""
	Create offspring of sequence1 and sequence2.

	This function creats and offspring of sequence1 and sequence2 by
	cutting them at a random point and merging the two ends.

	Parameters
	----------
	sequence1 : cosgen.sequence.Sequence
	    Parent one.
	sequence2 : cosgen.sequence.Sequence
	    Parent two.

	Returns
	-------
	cosgen.sequence.Sequence
	    Offspring for the two sequences given.
	"""
	length = min(len(sequence1.l),len(sequence2.l))
	r = random.randrange(length)
	l = np.append(sequence1.l[:r],sequence2.l[r:])
	nstimtypes = max(sequence1.nstimtypes,sequence2.nstimtypes)
	return sequence.Sequence(nstimtypes=nstimtypes,l=l)
