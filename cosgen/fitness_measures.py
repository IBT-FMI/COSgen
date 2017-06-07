import numpy as np

class OptimalityError(Exception):
	pass

def estimator_variance(sequence, model, optimality, contrast=None):

	X = model.design_matrix(sequence)
	
	covariance_beta = model.cov_beta(X)

	if contrast is not None:
		covariance_beta = np.dot(contrast,np.dot(covariance_beta,contrast.transpose()))
	if optimality=='a':
		return np.trace(covariance_beta)
	elif optimality=='d':
		return np.linalg.det(covariance_beta)
	else:
		raise OptimalityError('Unknow optimality {0}. Possible choices are "a","d"'.format(optimality))
