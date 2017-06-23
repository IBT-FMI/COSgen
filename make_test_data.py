#!/usr/bin/env python3

import nibabel as nib
import argh
import os.path
import numpy as np

def main(source,destination,sequence,hrf,voxels,scale=1000):
	seq = np.load(sequence)
	hrf = np.load(hrf)
	conv = np.convolve(seq,hrf)
	func_data = conv/conv.max()/scale
	roi = nib.load(voxels)
	vox = np.nonzero(roi.get_data())
	img = nib.load(os.path.expanduser(source))
	data = img.get_data()
	conv = np.convolve(seq,hrf)
	func_data = conv/conv.max()/scale
	func_data = func_data[0:img.shape[3]]
	data[vox] += func_data
	nib.save(img,os.path.expanduser(destination))

if __name__ == '__main__':
	argh.dispatch_command(main)
