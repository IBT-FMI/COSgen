import numpy as np
import os
import random

class BlockSizeError(Exception):
	pass

class Sequence:
	def __init__(self, seqlen=None, nstimtypes=1, seqtype='random', l=None, amplitudes=None, block_size=None, block_sigma=0, gap_size=None, gap_sigma=0):
		self.fitness = np.nan
		self.nstimtypes = nstimtypes
		if l is not None:
			self.l = np.array(l)
			self.seqlen = len(l)
			self.amplitudes = np.zeros(self.seqlen)
			idx = np.nonzero(self.l)
			if type(amplitudes)==str and amplitudes=='random':
				self.amplitudes[idx] = np.random.choice(np.linspace(0,1,256),len(idx))
			elif amplitudes is not None:
				self.amplitudes = np.array(amplitudes)
			else:
				self.amplitudes[idx]= np.ones(len(idx))
			if self.nstimtypes == 1:
				self.nstimtypes = len(set(l))-1
			if self.nstimtypes == 0:	#if list only containes 0 events assume nstimtypes is 1
				self.nstimtypes = 1
		elif seqtype=='random':
			self.l = np.random.randint(nstimtypes+1, size=seqlen)
			self.amplitudes = np.zeros(seqlen)
			idx = np.nonzero(self.l)
			if type(amplitudes)==str and amplitudes=='random':
				self.amplitudes[idx] = np.random.choice(np.linspace(0,1,256),len(idx))
			else:
				self.amplitudes[idx] = np.ones(len(idx))
			self.seqlen = seqlen
		elif seqtype=='block':
			if block_size is None:
				raise BlockSizeError("block_size must be set if seqtype is 'block'")
			if gap_size is None:
				gap_size = block_size
			self.l = np.empty(seqlen)
			self.amplitudes = np.zeros(seqlen)
			self.seqlen = seqlen
			position=0
			if type(amplitudes)==str and  amplitudes=='random':
				if type(block_size) != list:
					while position<seqlen:
						block_len = max(int(np.random.normal(block_size,block_sigma)),1)
						self.l[position:min(position+block_len,seqlen)] = random.randint(1,nstimtypes)
						self.amplitudes[position:min(position+block_len,seqlen)] = np.random.choice(np.linspace(0,1,256),1)
						block_gap = max(int(np.random.normal(gap_size,gap_sigma)),1)
						self.l[position+block_len:min(position+block_len+block_gap,seqlen)] = 0
						position += block_len+block_gap
				else:
					while position<seqlen:
						block_len = np.random.choice(block_size)
						self.l[position:min(position+block_len,seqlen)] = random.randint(1,nstimtypes)
						self.amplitudes[position:min(position+block_len,seqlen)] = np.random.choice(np.linspace(0,1,256),1)
						block_gap = max(int(np.random.normal(gap_size,gap_sigma)),1)
						self.l[position+block_len:min(position+block_len+block_gap,seqlen)] = 0
						position += block_len+block_gap
			else:
				if type(block_size) != list:
					while position<seqlen:
						block_len = max(int(np.random.normal(block_size,block_sigma)),1)
						self.l[position:min(position+block_len,seqlen)] = random.randint(1,nstimtypes)
						self.amplitudes[position:min(position+block_len,seqlen)] = 1
						block_gap = max(int(np.random.normal(gap_size,gap_sigma)),1)
						self.l[position+block_len:min(position+block_len+block_gap,seqlen)] = 0
						position += block_len+block_gap
				else:
					while position<seqlen:
						block_len = np.random.choice(block_size)
						self.l[position:min(position+block_len,seqlen)] = random.randint(1,nstimtypes)
						self.amplitudes[position:min(position+block_len,seqlen)] = 1
						block_gap = max(int(np.random.normal(gap_size,gap_sigma)),1)
						self.l[position+block_len:min(position+block_len+block_gap,seqlen)] = 0
						position += block_len+block_gap
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

	def add_baseline(self,length,position='initial'):
		if position=='initial':
			self.l = np.append(np.zeros(length),self.l)
			self.seqlen = len(self.l)
		elif position=='terminal':
			self.l = np.append(self.l,np.zeros(length))
			self.seqlen = len(slef.l)
		else:
			raise ValueError("baseline position can only be 'initial' or 'terminal'. {0} was given.".format(position))

	def dump(self, path, index=0, TR=1, name='sequence'):
		path = os.path.expanduser(path)
		np.save(os.path.join(path,name+str(index)+'.npy'),self.l)
		with open(os.path.join(path,name+str(index)+'.tsv'),'w+') as f:
			f.write('onset\tduration\tfrequency\tpulse_width\tamplitude\tout_channel\n')
			blocks = self.get_block_representation()
			for i in blocks:
				f.write(str(i[0]*TR)+'\t'+str((i[1]-i[0])*TR)+'\t20.0\t0.005\t'+str(self.amplitudes[i[0]])+'\t1\n')

def estimate_optimal_block_size(seqlen, fc):
	blockseqs = [Sequence(seqlen,seqtype='block',block_size=i) for i in range(1,seqlen+1)]
	blockfitnesses = [fc.evaluate_fitness(blockseqs[i]) for i in range(seqlen)]
	bestblocksize = np.argmax(blockfitnesses)+1
	return bestblocksize
