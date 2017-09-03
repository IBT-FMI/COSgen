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

		This class is the base class can be used to create
		customized models for the calculation of fitness
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
		Return design matrix for a given `sequence`.

		This method executes the `design_matrix_func` given in
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

	def cov_beta(self, X):
		"""
		Return covariance matrix for a given design matrix `X`.

		This method executes the `cov_beta_func` given in
		the initialisation of the object. The parameter types
		and return types depend on the particular function.

		Parameters
		----------
		X : design matrix
		    Design matrix for which the covariance matrix is to be
		    calculated.

		Returns
		-------
		covariance matrix
		    Covariance matrix for the given design matrix.
		"""
		return self.cov_beta_func(X)



class EstimationModel(Model):
	"""
	This class implements a model for estimating the Haemodynamic Response Function (HRF).

	The model employs pre-whitening to account for
	autocorrelation for the errors. Either `whitening_mat` or
	`err_cov_mat` must be given.

	Parameters
	----------
	basis_set : numpy array
	    Array with hrf basis vectors as rows.
	whitening_mat : numpy matrix, optional
	    Whitening matrix.
	err_cov_mat : numpy matrix, optional
	    Error covariance matrix.
	filterfunc : function
	    Filter function takes numpy array as input and returns filtered
	    numpy array (c.f. :func:`~cosgen.models.gaussian_highpass`)
	convolution_func : function
	    Function used for convolution of HRF and sequence. Can be
	    changed in order to correct for non-linearity.
	extra_evs : array-like object
	    Extra explanatory variables in form of a 2D array-like object
	    with regressors as columns. Shapes is
	    (number of extra EVs, sequence length). If None, a constant regressor
	    is used for baseline correction.
	"""

	def __init__(self, basis_set, whitening_mat=None, err_cov_mat=None, filterfunc=lambda x: x, convolution_func=np.convolve, extra_evs=None):
		self.basis_set = basis_set
		for i in range(len(basis_set)):
			self.basis_set[i]=self.basis_set[i]/float(self.basis_set[i].max())
		self.filterfunc = filterfunc
		self.convolution_func = convolution_func
		if whitening_mat is not None:
			self.whitening_mat = np.matrix(whitening_mat)
		elif err_cov_mat is not None:
			#L = np.linalg.cholesky(err_cov_mat)
			#self.whitening_mat = np.matrix(np.linalg.inv(L))
			self.whitening_mat = get_whitening_mat(err_cov_mat)
		else:
			raise AttributeError("Either 'whitening_mat or 'err_cov_mat' must be given.")
		if extra_evs is None:
			self.extra_evs = np.empty((len(self.whitening_mat),1))
			self.extra_evs[:,0] = np.ones(len(self.whitening_mat))
			#self.extra_evs = np.array([[]])
		else:
			self.extra_evs = extra_evs
		self.n_extra_evs = self.extra_evs.shape[1]

	def design_matrix(self, sequence):
		"""
		Calculate design matrix.

		This method calculates the design matrix for a given
		`sequence`. Columns of the design matrix are a constant
		(ones), a linearly increasing time course, and the convolution of the
		basis vectors with the `sequence`.

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
			idx = sequence.l == i
			tmp = np.zeros(sequence.seqlen)
			tmp[idx] = sequence.amplitudes[idx]
			for j in range(lb):
				conv = self.convolution_func(tmp, self.basis_set[j])
				DM[:, self.n_extra_evs + lb * (i-1) + j] = self.filterfunc(conv[0:ls])
		return DM

	def cov_beta(self, X):
		"""
		Calculate covariance of estimators (betas).

		This method calculates the covariance matrix of the
		estimators for a given design matrix `X`. It employs
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
		Z = np.dot(self.whitening_mat, X)
		try:
			Zpinv = np.linalg.pinv(Z)
		except np.linalg.linalg.LinAlgError:
			print('X:')
			print(X)
			print('K:')
			print(self.whitening_mat)
			print('Z:')
			print(Z)
		return np.dot(Zpinv, np.transpose(Zpinv))

class DetectionModel(Model):
	"""
	This class implements a model for detecting specific
	contrasts for a given/known Haemodynamic Response Function (HRF).

	Unless otherwise specified the parameters and returns are the same
	as for `EstimationModel`.

	Parameters
	----------
	hrf : np.array
	    HRF values at multiples of TR.
	"""

	def __init__(self, hrf, whitening_mat=None, err_cov_mat=None, filterfunc=lambda x: x, convolution_func=np.convolve, extra_evs=None):
		self.hrf = hrf/float(hrf.max())
		self.filterfunc = filterfunc
		self.convolution_func = convolution_func
		if whitening_mat is not None:
			self.whitening_mat = np.matrix(whitening_mat)
		elif err_cov_mat is not None:
			#L = np.linalg.cholesky(err_cov_mat)
			#self.whitening_mat = np.matrix(np.linalg.inv(L))
			self.whitening_mat = get_whitening_mat(err_cov_mat)
		else:
			raise AttributeError("Either 'whitening_mat or 'err_cov_mat' must be given.")
		if extra_evs is None:
			self.extra_evs = np.empty((len(self.whitening_mat),1))
			self.extra_evs[:,0] = np.ones(len(self.whitening_mat))
			#self.extra_evs = np.array([[]])
		else:
			self.extra_evs = extra_evs
		self.n_extra_evs = self.extra_evs.shape[1]

	def design_matrix(self, sequence):
		"""
		Calculate design matrix.

		This method calculates the design matrix for a given
		`sequence`. Columns of the design matrix are a constant
		(ones), a linearly increasing time course, and the convolution of `hrf`
		with the sequence.

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
		DM = np.empty((ls, self.n_extra_evs+sequence.nstimtypes))
		DM[:,0:self.n_extra_evs]=self.extra_evs
		for i in range(1, sequence.nstimtypes+1):
			idx = sequence.l == i
			tmp = np.zeros(sequence.seqlen)
			tmp[idx] = sequence.amplitudes[idx]
			try:
				conv = self.convolution_func(tmp, self.hrf)
			except TypeError:
				conv = self.convolution_func(tmp)
			DM[:,self.n_extra_evs + i-1 ] = orthogonalize(self.extra_evs, self.filterfunc(conv[0:ls]))
		return DM

	def cov_beta(self, X):
		"""
		Calculate covariance of estimators (betas).

		This method calculates the covariance matrix of the
		estimators for a given design matrix `X`. It employs
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
		Z = np.dot(self.whitening_mat, X)
		try:
			Zpinv = np.linalg.pinv(Z)
		except np.linalg.linalg.LinAlgError:
			print('X:')
			print(X)
			print('K:')
			print(self.whitening_mat)
			print('Z:')
			print(Z)
		return np.dot(Zpinv, np.transpose(Zpinv))


#def get_canonical_basis_set(TR, length, order):
#	#TODO implement properly
#	pass
#
#def get_gamma_basis_set(TR, length, order, a1, b1, a2, b2, c):
#	#TODO implement properly
#	t = range(0,length*TR,TR)
#	basis = []

def get_FIR_basis_set(length):
	"""
	Return Finite Response basis set.

	Parameters
	----------
	length : int
	    Length of Basis set.
	"""
	return np.identity(length)

#def get_bspline_basis_set(TR, length, order):
#	#TODO implement properly
#	pass
#
#def get_fourier_basis_set(TR, length, order):
#	#TODO implement properly
#	pass
#
#def get_ICA_basis_set(TR, length, order):
#	#TODO implement properly
#	pass

def get_gamma_hrf(TR, length, a1=6, a2=16, a3=1, a4=1, a5=1, a6=0):
	"""
	Return Haemodynamic Response Function (HRF) that is difference of two gamma function.

	This function returns an HRF array constructed using the following
	formula\:

	..  math::  h(t) = \\frac{(t-a_6)^{a_1/a_3-1}\\exp(-(t-a_6)/a_3)}{\\Gamma (a_1/a_3)a_3^{a_1/a_3}} - a_5\\frac{(t-a_6)^{a_2/a_4-1}\\exp(-(t-a_6)/a_4)}{\\Gamma (a_2/a_4)a_4^{a_2/a_4}}

	Defaults are the same as those used by SPM. The maximum is normalized to 1.

	Parameters
	----------
	TR : float
	    Repetition time.
	length : int
	    Length of hrf.
	a1 : float
	    Function parameter of hrf. Time to peak.
	a2 : float
	    Function parameter of hrf.
	a3 : float
	    Function parameter of hrf.
	a4 : float
	    Function parameter of hrf.
	a5 : float
	    Function parameter of hrf. Time to onset.
	a6 : float
	    Function parameter of hrf. Time shift.

	Returns
	-------
	np.array
	    Array with hrf values at multiples of TR.
	"""
	t = np.linspace(0, (length-1)*TR, length) - a6
	const1 = a1/a3
	const2 = a2/a4
	denom1 = scipy.special.gamma(const1)*a3**(const1)
	denom2 = scipy.special.gamma(const2)*a4**(const2)
	hrf = t**(const1-1) * np.exp(-t/a3)/denom1 - a5 * t**(const2-1) * np.exp(-t/a4)/denom2
	return hrf/sum(hrf)

def get_ar1_cov(dim, phi):
	"""
	Return covariance matrix for first order autoregressive model.

	Parameters
	----------
	dim : int
	    Dimensions.
	phi : float
	    Autoregressive parameter.

	Returns
	-------
	np.matrix
	    Covariance matrix.
	"""
	return np.matrix(np.fromfunction(lambda i, j: phi**np.abs(i-j), (dim, dim)))

def get_whitening_mat(cov_mat, whitening_type='zca', epsilon=1e-10):
	"""
	Calculate whitening matrix from covariance matrix.

	Parameters
	----------
	cov_mat : 2d array-like object
	    Covariance matrix.
	whitening_type : string, optional
	    Type of whitening. Can be 'zca' or 'cholesky'.
	epsilon : float
	    Constant to avoid zero division when inverting eigenvalues.

	Returns
	-------
	np.matrix
	    Whitening matrix.
	"""
	if whitening_type == 'zca':
		U,S,V = np.linalg.svd(cov_mat)
		whitening_mat = np.matrix(np.dot(U, np.dot(np.diag(1./np.sqrt(S + epsilon)), U.T)))
		return whitening_mat
	elif whitening_type == 'cholesky':
		L = np.linalg.cholesky(cov_mat)
		whitening_mat = np.matrix(np.linalg.inv(L))
		return whitening_mat
	else:
		raise ValueError('Unknown whitening type "{0}".')


def get_autocorr_whitening_mat(acf, epsilon=1e-10):
	"""
	Return ZCA whitening matrix for a given autocorrelation function.

	Parameters
	----------
	acf : np.array
	    Autocorrelation function.
	epsilon : float
	    Constant added to 0 eigenvalues to avoid division by zero.

	Returns
	-------
	np.matrix
	    Whitening matrix.
	"""
	U,S,V = np.linalg.svd(scipy.linalg.toeplitz(acf))
	whitening_mat = np.matrix(np.dot(U, np.dot(np.diag(1./np.sqrt(S + epsilon)), U.T)))
	return whitening_mat

def plot_design_matrix(mat):
	"""
	Show gray scale plot of design matrix.

	Parameters
	----------
	mat : np.matrix
	    Design matrix.
	"""
	a = int(mat.shape[0]/mat.shape[1])
	fig, ax = plt.subplots(1)
	plt.imshow(np.repeat(mat, a, axis=1), cmap='gray')
	plt.title('Design matrix')
	ax.set_ylabel('time course [TR]')
	ax.set_xlabel('Regressors')
	ax.set_xticklabels([])
	plt.show()

def orthogonalize(A, v):
	"""
	Return vector orthogonalized with respect to column vectors of
	matrix `A`.

	Beware `A` must contain already orthogonalized column vectors!!

	Parameters
	----------
	A : 2-d array-like object
	    Basis vectors.
	v : 1-d array-like object
	    Vector to be orthogonalized.

	Returns
	-------
	1-d array-like object
	    Orthogonalized vector.
	"""
	if A.shape[1] == 0:
		return v
	v = v - A[:,0].dot(v)/A[:,0].dot(A[:,0]) * A[:,0]
	if A.shape[1] > 1:
		return orthogonalize(A[:,1:], v)
	else:
		return v

def gaussian_highpass(data, sigma=225):
	"""
	Return filtered data.

	This function filters an 1-d array with a Gaussian
	high-pass filter.

	Parameters
	----------
	data : 1-d array-like object
	     Input data.
	sigma : int
	    Standard deviation of Gaussian kernel.
	"""
	return data - gaussian_filter1d(data, sigma)

def squashing_function(array, max=2):
	"""
	Reduce all values in `array` bigger than `max` to `max`.

	Parameters
	----------
	array : array like
	    Array of values to be squashes.
	max : float
	    Maximum value used for squashing.
	"""
	idx = np.where(array > max)
	array[idx] = max
	return array

def tukey_taper(data, m=15):
	"""
	Apply Tukey window of length `m` to `data`.

	Parameters
	----------
	data : array like
	    Data to be windowed.
	m : int
	    Length of Tukey window.
	"""
	t=np.arange(0, m, 1)
	w = 0.5*(1+np.cos(np.pi*t/m))
	result = np.zeros(len(data))
	result[:m]=w*data[:m]
	return result

def volterra_corrected_convolution(seq, kernel1, kernel0=None, kernel2=None, const=-0.3):
	"""
	Compute response to sequence using a Volterra series.

	Calculates the response to a sequence using a second order
	Volterra series.

	.. math:: y(t) = h_0 + \sum_{i=0}^Th_1(i)seq(t-i) + \sum_{i=0}^T\sum_{j=0}^Th_2(i,j)seq(t-i)seq(t-j)

	For a very simple non-linearity correction: kernel0=0,
	kernel1=hrf, kernel2=-hrf*hrf.T.

	Parameters
	----------
	seq : array like
	    Sequence.
	kernel1 : array like
	    First order Volterra kernel.
	kernel0 : float, optional
	    Zeroth order Volterra kernel.
	    Default is 0.
	kernel2 : array like 2D, optional
	    Second order Volterra kernel.
	    Default is cont*kernel1*kernel1.T.
	const : float, optional
	    Constant for construction of kernel2.
	    Default is -0.3.
	"""
	if kernel0 is None:
		kernel0=0
	if kernel2 is None:
		colv = np.reshape(kernel1,(len(kernel1),1))
		rowv = np.reshape(kernel1,(1,len(kernel1)))
		kernel2=np.dot(const*colv,rowv)
	normalresponse = np.convolve(seq, kernel1)
	conv2d = np.zeros(len(seq)+kernel2.shape[0]-1)
	for i in range(kernel2.shape[0]):
		tmp = np.convolve(kernel2[i,:],seq)
		for t in range(i,len(seq)+i):
			conv2d[t] += seq[t-i]*tmp[t]
	volterra_corrected = kernel0 + normalresponse + conv2d
	return volterra_corrected
