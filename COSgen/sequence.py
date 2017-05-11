import datetime
import random

class sequence:
	def __init__(self,seqlen,nstimtypes,seqtype='random'):
		if seqtype=='random':
			self.l = random.sample(range(nstimtypes),seqlen)
		elif seqtype=='block':
			pass
		elif seqtype=='m':
			pass
		self.fitness = None
	def dump(self,path):
		path = path + '/{:%Y%m%d%H%M%S}.tsv'.format(datetime.datetime.now())
		with open(path,'w+') as f:
			f.write('hello')
