"""This file contains the function necessary to generate immigrants.
Any alternative functions should be added here."""
try:
	import cosgen.sequence as sequence
except ImportError:
	import sequence

import numpy as np

def generate_immigrants(nimmigrants, seqlen, nstimtypes, block_size, gap_size, cross_over_fct, block_sigma=None, gap_sigma=0, amplitudes=None):
	"""
	Generate immigrants.

	This function generates 3 times `nimmigrants` sequences partially
	consisting of block and a random sequences as well as a mixture
	of both types.

	Parameters
	----------
	nimmigrants : int
	    Number of sequences to be generated.
	seqlen : int
	    Length of the sequences to be generated.
	nstimtypes : int
	    Number of stimulus types of the sequences.
	block_size : int
	    Mean size of the blocks in the block sequence part.
	gap_size : int
	    Mean size of the gaps in the block sequence part.
	cross_over_fct : function
	    Function taking two sequences as parameters and returning one.
	    (e.g. cosgen.cross_over.cross_over)
	block_sigma : float
	    Standard deviation of the block sizes.
	gap_sigma : float
	    Standard deviation of the gap sizes.
	amplitudes : string
	    'random' for randomly assigned stimulus amplitudes and None for fixed
	    amplitudes.

	Returns
	-------
	list of cosgen.sequence.Sequence
	    List of sequence according to parameters.
	"""
	immigrants = []
	immigrants.extend(generate_block_immigrants(nimmigrants, seqlen, nstimtypes, block_size, gap_size, block_sigma, gap_sigma, amplitudes))
	for i in range(nimmigrants):
		randseq = sequence.Sequence(seqlen, nstimtypes, seqtype='random',amplitudes=amplitudes)
		blockseq = sequence.Sequence(seqlen, nstimtypes, seqtype='block', block_size=block_size, gap_size=gap_size, block_sigma=block_sigma, gap_sigma=gap_sigma, amplitudes=amplitudes)
		if np.random.rand()<0.5:
			seq = cross_over_fct(randseq, blockseq)
		else:
			seq = cross_over_fct(blockseq, randseq)
		immigrants.append(seq)
		immigrants.append(sequence.Sequence(seqlen, nstimtypes, seqtype='random',amplitudes=amplitudes))
	return immigrants

def generate_block_immigrants(nimmigrants, seqlen, nstimtypes, block_size, gap_size, block_sigma=None, gap_sigma=0, amplitudes=None):
	"""
	Generate block type immigrants.

	This function generates `nimmigrants` block design sequences.

	Parameters
	----------
	nimmigrants : int
	    Number of sequences to be generated.
	seqlen : int
	    Length of the sequences to be generated.
	nstimtypes : int
	    Number of stimulus types of the sequences.
	block_size : int
	    Mean size of the blocks in the block sequence part.
	gap_size : int
	    Mean size of the gaps in the block sequence part.
	block_sigma : float
	    Standard deviation of the block sizes.
	gap_sigma : float
	    Standard deviation of the gap sizes.
	amplitudes : string
	    'random' for randomly assigned stimulus amplitudes and None for fixed
	    amplitudes.

	Returns
	-------
	list of cosgen.sequence.Sequence
	    List of sequence according to parameters.
	"""
	immigrants = []
	for i in range(nimmigrants):
		immigrants.append(sequence.Sequence(seqlen, nstimtypes, seqtype='block', amplitudes=amplitudes, block_size=block_size, block_sigma=block_sigma, gap_size=gap_size, gap_sigma=gap_sigma))
	return immigrants
