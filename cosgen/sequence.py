import datetime
import random
import os

class sequence:
	def __init__(self,seqlen,nstimtypes,seqtype='random'):
		if seqtype=='random':
			self.l = [ random.randrange(0,nstimtypes+1,1) for _ in range(seqlen) ]
		elif seqtype=='block':
			pass
		elif seqtype=='m':
			pass
		self.fitness = None
	def dump(self,path):
		os.makedirs(path, exist_ok=True)
		path = path + '/{:%Y%m%d%H%M%S}.tsv'.format(datetime.datetime.now())
		with open(path,'w+') as f:
			f.write('hello')
