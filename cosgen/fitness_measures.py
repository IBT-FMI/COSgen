"""This file contains fitness measure functions."""
import numpy as np

class OptimalityError(Exception):
	"""Error raised when a optimality is not 'a' or 'd'."""
	pass

def estimator_variance(sequence, model, optimality, contrast=None):
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
	    (determinat).

	Retruns
	-------
	float
	    optimality value
	"""
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
