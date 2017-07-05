import numpy as np
import scipy.linalg
import scipy.special
from scipy import signal
from scipy.ndimage.filters import gaussian_filter1d
try:
	from matplotlib import pyplot as plt
except:
	pass

class Model:

	def __init__(self, design_matrix_func, cov_beta_func):
		"""
		Base class for models.

		This class is the base class can be used to creat 
		coustomized models for the calculation of fitness 
		measures.

		Parameters
		----------
		design_matrix_func : function
		    Function for generation of design matrix from 
		    a sequence.
		cov_beta_func : function
		    Function returning the covariance matrix of the
		    estimators (betas) for a given sequence.
		"""
		self.design_matrix_func = design_matrix_func
		self.cov_beta_func = cov_beta_func

	def design_matrix(self, sequence):
		"""
		Retrun design matrix for a given sequence.

		This method execute the 'design_matrix_func' given in 
		the initialisation of the object. The parameter types
		and return types depend on the particular function.

		Parameters
		----------
		sequence
		    Sequence for which the design matrix is to be 
		    calculated.

		Returns
		-------
		design matrix
		    Design matrix for the given sequence.
		"""
		return self.design_matrix_func(sequence)

	def cov_beta(self,X):
		"""
		Retrun covarinace matrix for a given design matrix.

		This method execute the 'cov_beta_func' given in 
		the initialisation of the object. The parameter types
		and return types depend on the particular function.

		Parameters
		----------
		design matrix
		    Design matrix for which the covariance matrix is to be
		    calculated.

		Returns
		-------
		covarince matrix
		    Covarnace matrix for the given design matrix.
		"""
		return self.cov_beta_func(X)



class EstimationModel(Model):
	"""
	This class implements a model for estimating the hrf.

	The model employes pre-whitening to account for 
	autocorrelation for the errors. Either 'whitening_mat or 
	'err_cov_mat' must be given.

	Parameters
	----------
	basis_set : numpy array
	    Array with hrf basis vetors as rows.
	whitening_mat : numpy matrix, optional
	    Whitening matrix.
	err_cov_mat : numpy matrix, optional
	    Error covariance matrix.
	filterfunc : function
	    Filter function takes numpy array as input and returns filtered 
	    numpy array (c.f. :func:`~cosgen.models.gaussian_highpass`)
	extra_evs : array-like object
	    Extra explenatory variables in form of a 2D array-like object 
	    with regressors as collumns. Shapes is 
	    (number of extra evs, sequence length).
	"""

	def __init__(self, basis_set, whitening_mat=None, err_cov_mat=None, filterfunc=lambda x: x, extra_evs=None):
		self.basis_set = basis_set
		self.filterfunc = filterfunc
		if whitening_mat is not None:
			self.whitening_mat = whitening_mat
		elif err_cov_mat is not None:
			L = np.linalg.cholesky(err_cov_mat)
			self.whitening_mat = np.linalg.inv(L)
		else:
			raise AttributeError("Either 'whitening_mat or 'err_cov_mat' must be given.")
		if extra_evs is None:
			self.extra_evs = np.array([[]])
		else:
			self.extra_evs = extra_evs
		self.n_extra_evs = self.extra_evs.shape[1]

	def design_matrix(self, sequence):
		"""
		Calculate design matrix.

		This method calculates the desing matrix for a given
		sequence. Colums of the desing matrix are a constant 
		(ones) a linear time course and the convolution of the
		basis vetors with the sequence.

		Parameters
		----------
		sequence : cosgen.sequence.Sequence
		    Sequence for which the design matrix is calculated.

		Returns
		-------
		numpy matrix
		    Design matrix.
		"""
		lb = len(self.basis_set)
		ls = len(sequence.l)
		DM = np.empty(( ls, self.n_extra_evs + sequence.nstimtypes * lb ))
		DM[:,0:self.n_extra_evs]=self.extra_evs
		for i in range(1, sequence.nstimtypes+1):
			for j in range(lb):
				DM[:, self.n_extra_evs + lb * (i-1) + j] = self.filterfunc(np.convolve(sequence.l == i, self.basis_set[j])[0:ls])
		return DM

		DM = np.empty((ls,self.n_extra_evs+sequence.nstimtypes))
		DM[:,0:self.n_extra_evs]=self.extra_evs
		X = np.array([sequence.l == i for i in range(1,sequence.nstimtypes+1)],dtype=int)
		DM[:,self.n_extra_evs:] = np.transpose(np.apply_along_axis(lambda m: orthogonalize(self.extra_evs,self.filterfunc(np.convolve(m,self.hrf)[0:ls])), axis=1, arr=X))
		return DM

	def cov_beta(self, X):
		"""
		Calculate covariance of estimators (betas).

		This method calculated the covariance matrix of the 
		estimators for a given design matrix. It employs 
		pre-whitening.

		Parameters
		----------
		X : numpy matrix
		    Design matrix.

		Returns
		-------
		numpy matrix
		    Covariance matrix of beta.
		"""
		#This is only for pre-whitening and not precoloring
		Z = self.whitening_mat*X
		try:
			Zpinv = np.linalg.pinv(Z)
		except np.linalg.linalg.LinAlgError:
			print('X:')
			print(X)
			print('K:')
			print(self.whitening_mat)
			print('Z:')
			print(Z)
		return Zpinv * np.transpose(Zpinv)	

class DetectionModel(Model):

	def __init__(self, hrf, whitening_mat=None, err_cov_mat=None, filterfunc=lambda x: x, extra_evs=None):
		"""
		This class implements a model for detecting specific 
		constrasts for a given/known hrf.

		The model employes pre-whitening to account for 
		autocorrelation for the errors. Either 'whitening_mat or 
		'err_cov_mat' must be given.

		Parameters
		----------
		hrf : numpy array
		    Array with hrf values at multiples of TR.
		whitening_mat : numpy matrix, optional
		    Whitening matrix.
		err_cov_mat : numpy matrix, optional
		    Error covariance matrix.
		"""
		self.hrf = hrf
		self.filterfunc = filterfunc
		if whitening_mat is not None:
			self.whitening_mat = whitening_mat
		elif err_cov_mat is not None:
			L = np.linalg.cholesky(err_cov_mat)
			self.whitening_mat = np.linalg.inv(L)
		else:
			raise AttributeError("Either 'whitening_mat or 'err_cov_mat' must be given.")
		if extra_evs is None:
			self.extra_evs = np.array([[]])
		else:
			self.extra_evs = extra_evs
		self.n_extra_evs = self.extra_evs.shape[1]

	def design_matrix(self, sequence):
		"""
		Calculate design matrix.

		This method calculates the desing matrix for a given
		sequence. Colums of the desing matrix are a constant 
		(ones) a linear time course and the convolution of the
		hrf with the sequence.

		Parameters
		----------
		sequence : cosgen.sequence.Sequence
		    Sequence for which the design matrix is calculated.

		Returns
		-------
		numpy matrix
		    Design matrix.
		"""
		ls = len(sequence.l)
		DM = np.empty((ls,self.n_extra_evs+sequence.nstimtypes))
		DM[:,0:self.n_extra_evs]=self.extra_evs
		X = np.array([sequence.l == i for i in range(1,sequence.nstimtypes+1)],dtype=int)
		DM[:,self.n_extra_evs:] = np.transpose(np.apply_along_axis(lambda m: orthogonalize(self.extra_evs,self.filterfunc(np.convolve(m,self.hrf)[0:ls])), axis=1, arr=X))
		return DM

	def cov_beta(self, X):
		"""
		Calculate covariance of estimators (betas).

		This method calculated the covariance matrix of the 
		estimators for a given design matrix. It employs 
		pre-whitening.

		Parameters
		----------
		X : numpy matrix
		    Design matrix.

		Returns
		-------
		numpy matrix
		    Covariance matrix of beta.
		"""
		#This is only for pre-whitening and not precoloring
		Z = self.whitening_mat*X
		try:
			Zpinv = np.linalg.pinv(Z)
		except np.linalg.linalg.LinAlgError:
			print('X:')
			print(X)
			print('K:')
			print(self.whitening_mat)
			print('Z:')
			print(Z)
		return Zpinv * np.transpose(Zpinv)	


def get_canonical_basis_set(TR,length,order):
	pass

def get_gamma_basis_set(TR,length,order,a1,b1,a2,b2,c):
	t = range(0,length*TR,TR)
	basis = []

def get_FIR_basis_set(length):
	return np.identity(length)

def get_bspline_basis_set(TR,length,order):
	pass

def get_fourier_basis_set(TR,length,order):
	pass

def get_ICA_basis_set(TR,length,order):
	pass

def get_gamma_hrf(TR,length,a1,b1,a2,b2,c):
	t = np.linspace(0,(length-1)*TR,length)
	const1 = b1^(a1+1)
	const2 = b2^(a2+1)
	denom1 = scipy.special.gamma(a1+1)
	denom2 = c * scipy.special.gamma(a2+1)
	return const1 * t**a1 * np.exp(-b1*t)/denom1 - const2 * t**a2 * np.exp(-b2*t)/denom2

def get_ar1_cov(dim,phi):
	return np.matrix(np.fromfunction(lambda i, j: phi**np.abs(i-j), (dim, dim)))

def get_autocorr_whitening_mat(acf):
	return np.linalg.inv(np.linalg.cholesky(scipy.linalg.toeplitz(acf)))

def plot_design_matrix(mat):
	a = int(mat.shape[0]/mat.shape[1])
	fig, ax = plt.subplots(1)
	plt.imshow(np.repeat(mat,a,axis=1), cmap='gray')
	plt.title('Design matrix')
	ax.set_ylabel('time course [TR]')
	ax.set_xlabel('Regressors')
	ax.set_xticklabels([])
	plt.show()

def orthogonalize(A,v):
	"""A must contain already orthogonalized vector!!"""
	if A.shape[1] == 0:
		print('Not orthogonalized!')
		return v
	v = v - A[:,0].dot(v)/A[:,0].dot(A[:,0]) * A[:,0]
	if A.shape[1] > 1:
		return orthogonalize(A[:,1:],v)
	else:
		return v
	#coef = A.transpose().dot(v)/np.linalg.norm(A,axis=0)
	#return v - A.dot(coef)

def gaussian_highpass(data,sigma=225):
	return data - gaussian_filter1d(data,sigma)
