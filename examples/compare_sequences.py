#!/usr/bin/env python3

import custom_model
import cosgen
import numpy as np

block = cosgen.sequence.Sequence(1490,seqtype="block",block_size=196)
optseq = cosgen.sequence.Sequence(l=np.load('/home/wguest/.cosgen/sequence0.npy'))

model = cosgen.models.Model(custom_model.design_matrix, custom_model.cov_beta)
#dmblock = custom_model.design_matrix(block)
#dmoptseq = custom_model.design_matrix(optseq)

c = np.zeros(4)
c[0] = 0
c[1] = 0
c[2] = 1
c[3] = 0

print('Block fitness ', cosgen.fitness_measures.estimator_variance(block, model, 'd', np.matrix(c)))
print('Optimised seq fitness ', cosgen.fitness_measures.estimator_variance(optseq, model, 'd', np.matrix(c)))
