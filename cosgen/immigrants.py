"""This file contains the function neccessary to generate immigrants.
Any alternative functions should be added here."""
try:
	import cosgen.sequence as sequence
except ImportError:
	import sequence

def generate_immigrants(nimmigrants, seqlen, nstimtypes, block_size, cross_over_fct):
	"""
	Generate immigrants.

	This function generates 'nimmigrants' sequences partially consisting of a
	block and a random sequence.

	Parameters
	----------
	nimmigrants : int
	    Number of sequences to be generated.
	seqlen : int
	    Length of the sequences to be generated.
	nstimtypes : int
	    Number of stimulus types of the sequences.
	block_size : int
	    Size of the blocks in the block sequence part. Has to be a divisor
	    of the sequence length.
	cross_over_fct : function
	    Function taking two sequences as parameters and returning one.
	    (e.g. cosgen.cross_over.cross_over)

	Returns
	-------
	list of cosgen.sequence.Sequence
	    List of sequence according to parameters.
	"""
	imigrants = []
	for i in range(nimmigrants):
		randseq = sequence.Sequence(seqlen, nstimtypes, seqtype='random')
		blockseq = sequence.Sequence(seqlen, nstimtypes, seqtype='block', block_size=block_size)
		seq = cross_over_fct(randseq, blockseq)
		imigrants.append(seq)
	return imigrants
