try:
	import cosgen.sequence as sequence
except ImportError:
	import sequence

def generate_immigrants(nimigrants,seqlen,nstimtypes,cross_over_fct):
	imigrants = []
	for i in range(nimigrants):
		randseq = sequence.sequence(seqlen,nstimtypes,seqtype='random')
		blockseq = sequence.sequence(seqlen,nstimtypes,seqtype='block')
		seq = cross_over_fct(randseq,blockseq)
		imigrants.append(seq)
	return imigrants
