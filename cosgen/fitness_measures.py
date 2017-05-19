import numpy as np

class OptimalityError(Exception):
	pass

def estimator_variance(covariance_beta, optimality, contrast=None):
	if contrast is not None:
		covaraiance_beta = np.dot(contrast,np.dot(covariance_beta,contrast.transpose()))
	if optimality=='a':
		np.trace(covariance_beta)
	elif optimality=='d':
		np.linalg.det(covariance_beta)
	else:
		raise OptimalityError('Unknow optimality {0}. Possible choices are "a","d","c"'.format(optimality))
