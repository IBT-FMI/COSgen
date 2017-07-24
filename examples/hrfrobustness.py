from cosgen import models
from cosgen.sequence import estimate_optimal_block_size
from cosgen.function_crate import FunctionCrate
from cosgen.fitness_measures import estimator_variance
import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import functools


from scipy import stats, signal
def get_irf(
	a=2,
	b=10,
	resolution=1,
	):
	my_x = np.linspace(0,100,100*resolution)
	my_y = stats.beta.pdf(my_x/100, a, b)
	my_z = np.linspace(0,0,100*resolution)
	my_z[:20*resolution]=1

	irf = signal.deconvolve(my_y, my_z)[1]
	block_response = signal.convolve(irf,my_z)
	basis_function = signal.deconvolve(block_response, my_z)[1] #should be equal to irf

	return irf/irf.sum()
nx = 10
ny = 10
X = np.linspace(1,3,nx)
Y = np.linspace(9,11,ny)
hrfs = [[get_irf(x,y) for y in Y] for x in X]

#nx = 20
#ny = 5
#X = np.linspace(4,8,nx)
#Y = np.linspace(0,1,ny)
#hrfs = [[models.get_gamma_hrf(1,32,a1=x,a6=y) for y in Y] for x in X]

block_size = np.empty((nx,ny))

seql = 1490 #sequence length
ecm = np.identity(seql)	#error covarianve matrix
ecm = models.get_ar1_cov(seql,0.5)
extra_evs = np.empty((seql,2))
extra_evs[:,0]=np.ones(seql)	#constant confound regressor
extra_evs[:,1]=np.linspace(-0.5,0.5,seql)	#linear confound regressor
fc = FunctionCrate()
c = np.zeros(3)
c[0]=0
c[1]=0
c[2]=1
print(len(X))
print(len(Y))
print(len(hrfs))
print(len(hrfs[0]))
for i in range(nx):
	for j in range(ny):
		print('i',i)
		print('j',j)
		model = models.DetectionModel(hrfs[i][j],err_cov_mat=ecm, filterfunc=models.gaussian_highpass, extra_evs=extra_evs)
		fc.add_fitness_measure('est_var',functools.partial(estimator_variance,model=model,optimality='d',contrast=None))#np.matrix(c)))
		block_size[i,j]=estimate_optimal_block_size(seql,fc)

fig = plt.figure()
ax = fig.gca(projection='3d')
X, Y = np.meshgrid(X,Y)
print(X.shape)
print(Y.shape)
print(block_size.shape)
np.savez('figurearrays.npz',X,Y,block_size.transpose())
ax.plot_surface(X,Y,block_size.transpose())
plt.savefig('robustness.pdf')
