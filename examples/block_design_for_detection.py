#This code is similar to the code of the cosgen.sequence.estimate_optimal_block_size but includes plots

from cosgen.function_crate import FunctionCrate
from cosgen.fitness_measures import estimator_variance
from cosgen.sequence import Sequence
import cosgen.models

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.widgets import Slider
import functools
from scipy.ndimage.filters import gaussian_filter1d

seql = 1490

hrf = np.load('hrf.npy')

ecm = np.identity(seql)
ecm = cosgen.models.get_ar1_cov(seql,0.3)
extra_evs = np.empty((seql,2))
extra_evs[:,0]=np.ones(seql)	#constant confound regressor
extra_evs[:,1]=np.linspace(-0.5,0.5,seql)	#linear confound regressor
model = cosgen.models.DetectionModel(hrf,err_cov_mat=ecm, filterfunc=cosgen.models.gaussian_highpass, extra_evs=extra_evs)

blockseqs = [Sequence(seql,seqtype='block',block_size=i) for i in range(1,seql+1)]

fc = FunctionCrate()
c = np.zeros(3)
c[0]=-1
c[1]=-1
c[2]=1
fc.add_fitness_measure('est_var',functools.partial(estimator_variance,model=model,optimality='d',contrast=None))#np.matrix(c)))

blockfitnesses = [fc.evaluate_fitness(blockseqs[i]) for i in range(seql)]
bestblocksize = np.argmax(blockfitnesses)+1
print('Best block size: ',bestblocksize)
plt.plot(range(1,seql+1),blockfitnesses)
plt.xlabel("Block size in multiple of TR")
plt.ylabel("Efficiency in arbirtry units")
plt.title("Block size efficiency spectrum")
plt.show()
cosgen.models.plot_design_matrix(model.design_matrix(blockseqs[bestblocksize-1]))

index = bestblocksize-1
fig, ax = plt.subplots()
plt.xlabel("time [TR]")
plt.ylabel("Regressor [arbitrary units]")
plt.title('High-pass filtered and orthogonalized regressor')
plt.subplots_adjust(left=0.25, bottom=0.25)
t = range(seql)
amp = model.design_matrix(blockseqs[index])[:,2]
l, = plt.plot(t,amp)
plt.axis([0,seql,-1.01,1.01])

sliderax = plt.axes([0.25, 0.1, 0.65, 0.03], facecolor='lightgoldenrodyellow')

slider = Slider(sliderax, 'Block length', 0, seql-1, valinit=index)

def update(val):
	index = int(slider.val)
	dm = model.design_matrix(blockseqs[index])
	l.set_ydata(dm[:,2])
	fig.canvas.draw_idle()
slider.on_changed(update)

plt.show()
