class Model:

	def __init__(self, design_matrix_func, cov_beta_func):
		self.design_matrix_func = design_matrix_func
		self.cov_beta_func = cov_beta_func

	def design_matrix(self, sequence):
		return self.design_matrix_func(sequence)

	def cov_beta(self,X):
		return self.cov_beta_func(X)



class EstimationModel(Model):

	def __init__(self, basis_set, whitening_mat=None, err_cov_mat=None, noise_var=1, min_detectable_stim_dur=0):
		self.basis_set = basis_set
		self.noise_var = noise_var
		if whitening_mat is not None:
			self.whitening_mat = whitening_mat
		if err_cov_mat is not None:
			L = np.linalg.cholesky(err_cov_mat)
			try:
				self.whitening_mat = np.linalg.inv(L)
			except LinAlgError:
				self.whitening_mat = np.linalg.pinv(L)

	def design_matrix(self, sequence):
		pass

	def cov_beta(self, X):
		pass

class DetectionModel(Model):

	def __init__(self, hrf, whitening_mat=None, err_cov_mat=None, min_detectable_stim_dur=0):
		self.hrf = hrf
		if whitening_mat is not None:
			self.whitening_mat = whitening_mat
		if err_cov_mat is not None:
			L = np.linalg.cholesky(err_cov_mat)
			try:
				self.whitening_mat = np.linalg.inv(L)
			except LinAlgError:
				self.whitening_mat = np.linalg.pinv(L)

	def design_matrix(self, sequence):
		X = [np.array(sequence.l) == i for i in range(sequence.nstimtypes+1)]
		Z = np.transpose(np.apply_along_axis(lambda m: np.convolve(m,self.hrf), axis=1, arr=np.array(X).astype(int)))
		return np.c_[np.ones(len(Z)),np.c_[range(len(Z)),Z]] #add base line and linear trend/drift

	def cov_beta(self, X):
		pass


def get_canonical_basis_set(TR,length,order):
	pass

def get_gamma_basis_set(TR,length,order):
	pass

def get_FIR_basis_set(TR,length,order):
	pass

def get_bspline_basis_set(TR,length,order):
	pass

def get_fourier_basis_set(TR,length,order):
	pass

def get_ICA_basis_set(TR,length,order):
	pass
