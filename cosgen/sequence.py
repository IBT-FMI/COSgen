import numpy as np
import os
import random

class BlockSizeError(Exception):
	pass

class Sequence:
	def __init__(self, seqlen=None, nstimtypes=1, seqtype='random', l=None, block_size=None):
		self.fitness = None
		self.nstimtypes = nstimtypes
		if l is not None:
			self.l = np.array(l)
			self.seqlen = len(l)
			if self.nstimtypes == 1:
				self.nstimtypes = len(set(l))-1
			if self.nstimtypes == 0:	#if list only containes 0 events assume nstimtypes is 1
				self.nstimtypes = 1
		elif seqtype=='random':
			self.l = np.random.randint(nstimtypes+1,size=seqlen)
			self.seqlen = seqlen
		elif seqtype=='block':
			if block_size < 0 or seqlen % block_size != 0:
				raise BlockSizeError('block_size must be a positive divisor of seqlen. block_size={0} seqlen={1}'.format(block_size,seqlen))
			self.l = np.empty(seqlen)
			position=0
			while position<seqlen:
				self.l[position:position+block_size] = random.randint(0,nstimtypes)
				position += block_size
		elif seqtype=='m':
			pass
	def dump(self, path):
		np.save(os.path.join(path,'sequence.npy'),self.l)
		with open(os.path.join(path,'sequence.txt'),'a+') as f:
			f.write('Sequence: '+str(self.l)+'\n')
			f.write('Fitness: '+str(self.fitness)+'\n')

def estimate_optimal_block_size(hrf):
	#returns an estimate for the optimal value of a block design in multiples of TR (if hrf has not enough sample points this might not work properly?)
	ft_hrf = np.absolute(np.fft.rfft(hrf))
	ft_freq = np.fft.rfftfreq(hrf.size)
	idx = np.argmax(ft_hrf)
	return 1./ft_freq[idx]
