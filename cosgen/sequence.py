import numpy as np
import os
import random

class BlockSizeError(Exception):
	pass

class Sequence:
	def __init__(self, seqlen=None, nstimtypes=1, seqtype='random', l=None, amplitudes=None, block_size=None, block_sigma=0, gap_size=None, gap_sigma=0):
		"""
		Class for sequences.

		This class represents sequences in the optimization.

		Parameters
		----------
		seqlen : int, opt
		    Length of sequences.
		nstimtypes : int, opt
		    Number of stimulus types in a sequence.
		    Default is 1.
		seqtype : string, opt
		    Sequence type can be 'random' or 'block'.
		    Default is 'random'.
		l : list, opt
		    List of with representation of a sequence.
		    Can be used to transform [0,0,0,1,1,1,0,0,...] into an
		    object of sequence class.
		amplitudes : list, opt
		    Can be used to specify the amplitudes of `l`.
		block_size : int, opt
		    Mean block size of 'block' type sequence.
		block_sigma : float, opt
		    Standard deviation of block sizes in a 'block' style sequence.
		gap_size : int, opt
		    Mean gap size of 'block' type sequence.
		gap_sigma : float, opt
		    Standard deviation of gap sizes in a 'block' style sequence.
		"""
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
		"""
		Find blocks in a sequence represented as list.

		Returns
		-------
		numpy array
		    First column are start positions and second column are
		    end positions of blocks.
		"""
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
		"""
		Add a baseline to the beginning or end of the sequence.

		This function add `length` rest events to a sequence.

		Parameters
		----------
		length : int
		    Number of rest events.
		position : string
		    Can be 'initial' for adding the rest events at the
		    beginning or 'terminal' for adding them at the end.
		"""
		if position=='initial':
			self.l = np.append(np.zeros(length),self.l)
			self.seqlen = len(self.l)
		elif position=='terminal':
			self.l = np.append(self.l,np.zeros(length))
			self.seqlen = len(slef.l)
		else:
			raise ValueError("baseline position can only be 'initial' or 'terminal'. {0} was given.".format(position))

	def dump(self, path, index=0, TR=1, name='sequence'):
		"""
		Save sequence to file.

		This function save the sequence as a 'tsv' file and as a
		list in a '.npy' file.

		Parameters
		----------
		path : string
		    Path to folder where files are stored.
		index : int, opt
		    Number added to file name.
		TR : float, opt
		    Repetition time in seconds.
		    Default is 1s.
		name : string, opt
		    File name.
		    Default is 'sequence'.
		"""
		path = os.path.expanduser(path)
		np.save(os.path.join(path,name+str(index)+'.npy'),self.l)
		with open(os.path.join(path,name+str(index)+'.tsv'),'w+') as f:
			f.write('onset\tduration\tfrequency\tpulse_width\tamplitude\tout_channel\n')
			blocks = self.get_block_representation()
			for i in blocks:
				f.write(str(i[0]*TR)+'\t'+str((i[1]-i[0])*TR)+'\t20.0\t0.005\t'+str(self.amplitudes[i[0]])+'\t1\n')

def estimate_optimal_block_size(seqlen, fc):
	"""
	Find optimal block size.

	This function finds the optimal block size for a sequence of
	length `seqlen` using the fitness measures in the
	`FunctionCrate` object `fc`.

	Parameters
	----------
	seqlen : int
	    Length of the sequence.
	fc : FunctionCrate object
	    Object containing fitness measures.
	"""
	blockseqs = [Sequence(seqlen,seqtype='block',block_size=i) for i in range(1,seqlen+1)]
	blockfitnesses = [fc.evaluate_fitness(blockseqs[i]) for i in range(seqlen)]
	bestblocksize = np.argmax(blockfitnesses)+1
	return bestblocksize
