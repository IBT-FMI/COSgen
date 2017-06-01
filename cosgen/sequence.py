import datetime
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
			self.l = l
			self.seqlen = len(l)
			self.nstimtypes = len(set(l))-1
		elif seqtype=='random':
			self.l = np.random.randint(nstimtypes+1,size=seqlen)
			self.seqlen = seqlen
		elif seqtype=='block':
			if block_size < 0 or seqlen % block_size != 0:
				raise BlockSizeError('block_size must be a positive devisor of seqlen. block_size={0} seqlen={1}'.format(block_size,seqlen))
			self.l = np.empty(seqlen)
			position=0
			while position<seqlen:
				self.l[position:position+block_size] = random.randint(0,nstimtypes)
				position += block_size
		elif seqtype=='m':
			pass
	def dump(self, path):
		os.makedirs(path, exist_ok=True)
		path = path + '/{:%Y%m%d%H%M%S}.tsv'.format(datetime.datetime.now())
		with open(path,'w+') as f:
			f.write('hello')
