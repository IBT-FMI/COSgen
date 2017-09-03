"""This file contains fitness measure functions."""
import numpy as np

class OptimalityError(Exception):
	"""Error raised when a optimality is not 'a' or 'd'."""
	pass

def estimator_variance(sequence, model, optimality, contrast=None, normalization=1):
	"""
	The optimality of the estimator variances.

	This function calculates the a- or d-optimality value of the
	covariance matrix of the estimators.

	Parameters
	----------
	sequence : cosgen.sequence.Sequence
	    Sequence for which the covariance matrix is calculated.
	model : cosgen.models.Model
	    Model class providing functions for the construction of the
	    design matrix and covariance matrix.
	optimality : string
	    Can be 'a' for a-optimality (trace) or 'd' for d-optimality
	    (determinant).
	contrast : numpy matrix
	    Matrix containing contrast vectors as rows.
	Returns
	-------
	float
	    optimality value
	"""
	X = model.design_matrix(sequence)

	covariance_beta = model.cov_beta(X)

	if contrast is not None:
		contrast = np.matrix(contrast)
		covariance_beta = contrast*covariance_beta*contrast.transpose()
	if optimality=='a':
		return 1./np.trace(covariance_beta)/normalization
	elif optimality=='d':
		return 1./np.linalg.det(covariance_beta)/normalization
	else:
		raise OptimalityError('Unknown optimality {0}. Possible choices are "a","d"'.format(optimality))

def jitter(sequence, normalization=1):
	sp = abs(np.fft.fft(sequence.l))
	return sum(sp)/sp.max()/len(sp)/normalization
