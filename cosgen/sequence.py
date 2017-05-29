import datetime
import random
import os

class sequence:
	def __init__(self, seqlen=None, nstimtypes=1, seqtype='random', l=None):
		if l is not None:
			self.l = l
			self.seqlen = len(l)
		elif seqtype=='random':
			self.l = [ random.randrange(0,nstimtypes+1,1) for _ in range(seqlen) ]
			self.seqlen = seqlen
		elif seqtype=='block':
			pass
		elif seqtype=='m':
			pass
		self.fitness = None
		self.nstimtypes = nstimtypes
	def dump(self, path):
		os.makedirs(path, exist_ok=True)
		path = path + '/{:%Y%m%d%H%M%S}.tsv'.format(datetime.datetime.now())
		with open(path,'w+') as f:
			f.write('hello')
