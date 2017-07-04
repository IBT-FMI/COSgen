import numpy as np
import os
import random

class BlockSizeError(Exception):
	pass

class Sequence:
	def __init__(self, seqlen=None, nstimtypes=1, seqtype='random', l=None, block_size=None):
		self.fitness = np.nan
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
			if block_size is None:
				raise BlockSizeError("block_size must be set if seqtype is 'block'")
			#if block_size < 0 or seqlen % block_size != 0:
			#	raise BlockSizeError('block_size must be a positive divisor of seqlen. block_size={0} seqlen={1}'.format(block_size,seqlen))
			self.l = np.empty(seqlen)
			self.seqlen = seqlen
			position=0
			while position<seqlen:
				self.l[position:min(position+block_size,seqlen)] = random.randint(1,nstimtypes)
				self.l[position+block_size:min(position+2*block_size,seqlen)] = 0
				position += 2*block_size
		elif seqtype=='m':
			#TODO implement m sequences
			pass

	def get_block_representation(self):
		result = []
		start = -1
		for i in range(self.seqlen):
			if self.l[i] != 0 and start<0:
				start = i
			elif self.l[i] == 0 and start>=0:
				result.append([start, i])
				start = -1
		if self.l[self.seqlen-1] != 0:
			result.append([start, self.seqlen])
		return np.array(result)
		
	def dump(self, path, index=0, TR=1):
		path = os.path.expanduser(path)
		np.save(os.path.join(path,'sequence'+str(index)+'.npy'),self.l)
		with open(os.path.join(path,'sequence'+str(index)+'.tsv'),'w+') as f:
			f.write('onset\tduration\tstimulation_frequency\n')
			blocks = self.get_block_representation()
			for i in blocks:
				f.write(str(i[0]*TR)+'\t'+str((i[1]-i[0])*TR)+'\t20.0\n')

def estimate_optimal_block_size(seqlen, fc):
	blockseqs = [Sequence(seqlen,seqtype='block',block_size=i) for i in range(1,seqlen+1)]	
	blockfitnesses = [fc.evaluate_fitness(blockseqs[i]) for i in range(seqlen)]
	bestblocksize = np.argmax(blockfitnesses)+1
	return bestblocksize
