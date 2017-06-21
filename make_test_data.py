#!/usr/bin/env python3

import nibabel as nib
import argh
import os.path
import numpy as np

def main(source,destination,sequence,hrf,voxels,scale=1):
	seq = np.load(sequence)
	hrf = np.load(hrf)
	conv = np.convolve(seq,hrf)
	func_data = conv/conv.max()/scale
	vox = np.load(voxels)
	img = nib.load(os.path.expanduser(source))
	data = img.get_data()
	conv = np.convolve(seq,hrf)
	func_data = conv/conv.max()/scale
	func_data = func_data[0:img.shape[3]]
	for i in len(vox):
		data[vox[i,0],vox[i,1],vox[i,2]]+=func_data
	nib.save(img,os.path.expanduser(destination))

if __name__ == '__main__':
	argh.dispatch_command(main)
