import numpy as np
import scipy.special
from matplotlib import pyplot as plt

class Model:

	def __init__(self, design_matrix_func, cov_beta_func):
		self.design_matrix_func = design_matrix_func
		self.cov_beta_func = cov_beta_func

	def design_matrix(self, sequence):
		return self.design_matrix_func(sequence)

	def cov_beta(self,X):
		return self.cov_beta_func(X)



class EstimationModel(Model):

	def __init__(self, basis_set, whitening_mat=None, err_cov_mat=None):
		self.basis_set = basis_set
		if whitening_mat is not None:
			self.whitening_mat = whitening_mat
		elif err_cov_mat is not None:
			L = np.linalg.cholesky(err_cov_mat)
			self.whitening_mat = np.linalg.inv(L)

	def design_matrix(self, sequence):
		lb = len(self.basis_set)
		Xconv = np.empty(( len(sequence.l) + len(self.basis_set[0]) - 1, sequence.nstimtypes * lb ))
		Xconv = Xconv/Xconv.max()
		for i in range(1, sequence.nstimtypes+1):
			for j in range(lb):
				Xconv[:, lb * (i-1) + j] = np.convolve(sequence.l == i, self.basis_set[j])
		return np.matrix(np.c_[np.ones(len(Xconv)),np.c_[np.linspace(0,1,len(Xconv)),Xconv]]) #add base line and linear trend/drift

	def cov_beta(self, X):
		#This is only for pre-whitening and not precoloring
		Z = self.whitening_mat*X
		Zpinv = np.linalg.pinv(Z)
		return Zpinv * np.transpose(Zpinv)	

class DetectionModel(Model):

	def __init__(self, hrf, whitening_mat=None, err_cov_mat=None):
		self.hrf = hrf
		if whitening_mat is not None:
			self.whitening_mat = whitening_mat
		elif err_cov_mat is not None:
			L = np.linalg.cholesky(err_cov_mat)
			self.whitening_mat = np.linalg.inv(L)

	def design_matrix(self, sequence):
		X = np.array([sequence.l == i for i in range(1,sequence.nstimtypes+1)],dtype=int)
		Xconv = np.transpose(np.apply_along_axis(lambda m: np.convolve(m,self.hrf), axis=1, arr=X))
		Xconvmax = Xconv.max()
		if Xconvmax != 0:
			Xconv = Xconv/Xconvmax
		return np.matrix(np.c_[np.ones(len(Xconv)),np.c_[np.linspace(0,1,len(Xconv)),Xconv]]) #add base line and linear trend/drift

	def cov_beta(self, X):
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

def get_gamma_basis_set(TR,length,order,a1,b1,a2,b2i,c):
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

def plot_design_matrix(mat):
	plt.imshow(mat, cmap='gray')
	plt.show()
