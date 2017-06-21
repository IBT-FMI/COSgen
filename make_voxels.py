#!/usr/bin/env python3

import numpy as np

vox = np.zeros((1000,3))

idx = 0
for i in range(25,35):
	for j in range(45,55):
		for k in range(20,30):
			vox[idx,:] = [i,j,k]
			idx += 1
np.save('voxels.npy',vox)
