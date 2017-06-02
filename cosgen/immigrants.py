try:
	import cosgen.sequence as sequence
except ImportError:
	import sequence

def generate_immigrants(nimigrants, seqlen, nstimtypes, block_size, cross_over_fct):
	imigrants = []
	for i in range(nimigrants):
		randseq = sequence.Sequence(seqlen, nstimtypes, seqtype='random')
		blockseq = sequence.Sequence(seqlen, nstimtypes, seqtype='block', block_size=block_size)
		seq = cross_over_fct(randseq, blockseq)
		imigrants.append(seq)
	return imigrants
