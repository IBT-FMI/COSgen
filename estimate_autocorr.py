#!/usr/bin/env python3

import nibabel as nib
import numpy as np
import argh
import os.path
import math


def main(data="~/test_ni_data/ofM.dr/preprocessing/as_composite/sub-5703/ses-ofM/func/sub-5703_ses-ofM_trial-EPI_CBV_chr_longSOA.nii.gz", 
	voxels="~/ni_data/templates/DSURQEc_200micron_mask.nii.gz"):

	roi = nib.load(os.path.expanduser(voxels))
	vox = np.nonzero(roi.get_data())
	img = nib.load(os.path.expanduser(data))
	data = img.get_data()[vox]
	print(data.shape)
	a = np.empty((2,data.shape[1]-1))
	corrcoeff = 0
	n = data.shape[0]
	for i in range(data.shape[0]):
		a[0,:] = data[i,:data.shape[1]-1]
		a[1,:] = data[i,1:]
		tmp = np.corrcoef(a)[0,1]
		if math.isnan(tmp):
			#print(a)
			#print("Skipped the voxel shown above because it's temporal correlation is 'nan'.")
			n -= 1
			continue
		corrcoeff += tmp
	print("{0} of {1} voxels were used for the calculation.".format(n,data.shape[0]))
	print("The first order auto correlation parameter is approximately " + str(corrcoeff/n))

if __name__ == '__main__':
	argh.dispatch_command(main)
